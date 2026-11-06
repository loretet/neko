! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Defines wall model flux bc, extending Neumann.
module wall_model_flux
  use num_types, only : rp
  use coefs, only : coef_t
  use bc, only : bc_t, BC_NEUMANN
  use field, only : field_t
  use field_list, only : field_list_t
  use vector, only : vector_t
  use utils, only : neko_error, nonlinear_index
  use json_module, only : json_file
  use json_utils, only : json_get
  use math, only : masked_gather_copy_0
  use device_math, only : device_masked_gather_copy_0
  use device_neumann, only : device_neumann_apply_scalar
  use neko_config, only : NEKO_BCKND_DEVICE
  use, intrinsic :: iso_c_binding, only : c_ptr
  use time_state, only : time_state_t
  implicit none
  private

  !> Takes the output from the wall model when this runs for heat or scalars.
  !! A separate bc type is necessary because wall_model_bc is vector-only.
  !! If that is used for a scalar (which raises an error right now), 
  !! partial_init_base would look up scheme_name // "_mu" and "_rho", which
  !! do not exist for temperature/scalars. Even if those were fixed, you’d 
  !! get a 2nd MOST instance: double cost, double log, and "tau" overwritten.
  !! The type stores a separate dummy field `field_bc`, which is passed
  !! to the wall model and can be populated with its output flux. The
  !! boundary condition then gathers these values at the bc mask locations
  !! into a compact flux vector and applies the weak neumann contribution.
  type, public, extends(neumann_t) :: wall_model_flux_t
     !> A dummy field which can be manipulated by the user to set flux values.
     type(field_t) :: field_bc
     !> A field list, which just stores `field_bc`, for convenience.
     type(field_list_t) :: field_list
     !> Flux values gathered on boundary-mask ordering.
     type(vector_t) :: flux
     !> Function pointer to the user routine performing the update of the flux.
     procedure(wall_model_flux_update), nopass, pointer :: update => null()
   contains
     !> Constructor.
     procedure, pass(this) :: init => wall_model_flux_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          wall_model_flux_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => wall_model_flux_free
     !> Finalize.
     procedure, pass(this) :: finalize => wall_model_flux_finalize
     !> Apply scalar by adding weak neumann contribution.
     procedure, pass(this) :: apply_scalar => wall_model_flux_apply_scalar
     !> (No-op) Apply vector.
     procedure, pass(this) :: apply_vector => wall_model_flux_apply_vector
     !> (No-op) Apply vector (device).
     procedure, pass(this) :: apply_vector_dev => wall_model_flux_apply_vector_dev
     !> Apply scalar (device).
     procedure, pass(this) :: apply_scalar_dev => wall_model_flux_apply_scalar_dev
     !> Gather flux values at masked points.
     procedure, pass(this), private :: gather_flux => wall_model_flux_gather_flux

  end type wall_model_flux_t

  !> Abstract interface defining a neumann condition on a list of fields.
  abstract interface
     subroutine wall_model_flux_update(fields, bc, time)
       import field_list_t, wall_model_flux_t, time_state_t
       type(field_list_t), intent(inout) :: fields
       type(wall_model_flux_t), intent(in) :: bc
       type(time_state_t), intent(in) :: time
     end subroutine wall_model_flux_update
  end interface

  public :: wall_model_flux_update

contains
  !> Constructor.
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine wall_model_flux_init(this, coef, json)
    class(wall_model_flux_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: scalar_name

    call json_get(json, "scalar_field", scalar_name)
    call this%init_from_components(coef, scalar_name)
    if (allocated(scalar_name)) then
       deallocate(scalar_name)
    end if

  end subroutine wall_model_flux_init

  !> Constructor from components.
  !! @param[in] coef The SEM coefficients.
  !! @param[in] scalar_name Name of the temperature field.
  subroutine wall_model_flux_init_from_components(this, coef, scalar_name)
    class(wall_model_flux_t), intent(inout), target :: this
    type(coef_t), intent(in) :: coef
    character(len=*), intent(in) :: scalar_name

    call this%init_base(coef)
    this%bc_type = BC_NEUMANN

    call this%field_bc%init(this%dof, scalar_name)
    call this%field_list%init(1)
    call this%field_list%assign_to_field(1, this%field_bc)

  end subroutine wall_model_flux_init_from_components

  !> Destructor.
  subroutine wall_model_flux_free(this)
    class(wall_model_flux_t), target, intent(inout) :: this

    call this%free_base
    call this%field_bc%free()
    call this%field_list%free()
    call this%flux%free()

    if (associated(this%update)) then
       this%update => null()
    end if

  end subroutine wall_model_flux_free

  !> Gather field-defined values into compact boundary flux storage.
  subroutine wall_model_flux_gather_flux(this)
    class(wall_model_flux_t), intent(inout) :: this

    if (this%msk(0) .gt. 0) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_masked_gather_copy_0(this%flux(1), this%field_bc%x_d, &
               this%msk_d, this%field_bc%dof%size(), this%msk(0))
       else
          call masked_gather_copy_0(this%flux(1), this%field_bc%x, &
               this%msk, this%field_bc%dof%size(), this%msk(0))
       end if
    end if

  end subroutine wall_model_flux_gather_flux

  !> Apply scalar by adding weak neumann contribution.
  !! @param x Field to which the weak neumann contribution is added.
  !! @param n Size of the array `x`.
  !! @param time The current time state.
  subroutine wall_model_flux_apply_scalar(this, x, n, time, strong)
    class(wall_model_flux_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    integer :: i, m, k, facet
    integer :: idx(4)
    real(kind=rp) :: area
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (.not. strong_) then
      cp => neko_registry%get_field("temperature_cp")

       !$omp single
       if (.not. this%updated) then
          call this%update(this%field_list, this, time)
          call this%gather_flux()
          this%updated = .true.
       end if
       !$omp end single

       m = this%facet_node_msk(0)
       !$omp do
       do i = 1, m
          k = this%facet_node_msk(i)
          facet = this%facet(i)
          idx = nonlinear_index(k, this%coef%Xh%lx, this%coef%Xh%lx, &
               this%coef%Xh%lx)
          area = 0.0_rp
          select case (facet)
          case (1,2)
             area = this%coef%area(idx(2), idx(3), facet, idx(4))
          case (3,4)
             area = this%coef%area(idx(1), idx(3), facet, idx(4))
          case (5,6)
             area = this%coef%area(idx(1), idx(2), facet, idx(4))
          end select
          !$omp atomic
          x(k) = x(k) + this%flux%x(i) * area * cp * this%rho%x(i)
       end do
       !$omp end do
    end if

  end subroutine wall_model_flux_apply_scalar

  !> Apply scalar (device).
  !! @param x_d Device pointer to the field to update.
  !! @param time The current time state.
  !! @param strm Device stream.
  subroutine wall_model_flux_apply_scalar_dev(this, x_d, time, strong, strm)
    class(wall_model_flux_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (.not. strong_) then
       !$omp single
       if (.not. this%updated) then
          call this%update(this%field_list, this, time)
          call this%gather_flux()
          this%updated = .true.
       end if
       !$omp end single

       if (this%facet_node_msk(0) .gt. 0) then
          call device_neumann_apply_scalar(this%facet_node_msk_d, &
               this%facet_d, x_d, &
               this%flux%x_d, this%coef%area_d, this%coef%Xh%lx, &
               size(this%facet_node_msk), strm)
       end if
    end if

  end subroutine wall_model_flux_apply_scalar_dev

  !> (No-op) Apply vector.
  subroutine wall_model_flux_apply_vector(this, x, y, z, n, time, strong)
    class(wall_model_flux_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong

    call neko_error("wall_model_flux cannot apply vector BCs.")

  end subroutine wall_model_flux_apply_vector

  !> (No-op) Apply vector (device).
  subroutine wall_model_flux_apply_vector_dev(this, x_d, y_d, z_d, time, &
       strong, strm)
    class(wall_model_flux_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(c_ptr), intent(inout) :: y_d
    type(c_ptr), intent(inout) :: z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm

    call neko_error("wall_model_flux cannot apply vector BCs.")

  end subroutine wall_model_flux_apply_vector_dev

  !> Finalize.
  subroutine wall_model_flux_finalize(this)
    class(wall_model_flux_t), target, intent(inout) :: this

    call this%finalize_base()
    call this%flux%init(this%msk(0))

  end subroutine wall_model_flux_finalize

end module wall_model_flux
