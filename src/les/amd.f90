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
!
!> Implements `amd_t`, the anisotropic minimum-dissipation (AMD) LES model.
module amd
  use num_types, only : rp
  use les_model, only : les_model_t
  use field, only : field_t
  use fluid_scheme_base, only : fluid_scheme_base_t
  use json_utils, only : json_get_or_default, json_get, json_get_or_lookup, &
       json_get_or_lookup_or_default
  use json_module, only : json_file
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_memcpy, HOST_TO_DEVICE
  use utils, only : neko_error
  use map_1d, only : map_1d_t
  use amd_cpu, only : amd_compute_cpu
  use amd_device, only : amd_compute_device
  use registry, only : neko_registry
  use logger, only : LOG_SIZE, neko_log
  implicit none
  private

  !> Implements the anisotropic minimum-dissipation (AMD) LES model.
  !! The eddy viscosity and (optionally) the eddy diffusivity of a scalar are
  !! computed following Rozema et al. (2015) and Abkar & Moin (2017), Eqs.
  !! (17)-(18), including the buoyancy contribution for stratified flows.
  !! The scaled gradients use the filter-width tensor G of
  !! `les_model_compute_delta_tensor`, so that for any element shape
  !! (dk^ f)(dk^ g) = c * df/dx_a * G_ab * dg/dx_b.
  !! @note Reference DOIs: 10.1063/1.4928700, 10.1007/s10546-017-0288-4
  type, public, extends(les_model_t) :: amd_t
     !> Modified Poincare constant, defaults to 1/12.
     real(kind=rp) :: c
     !> Whether to use anisotropic filter widths (the filter-width tensor).
     logical :: if_aniso = .true.
     !> Whether to compute the eddy diffusivity of a scalar.
     logical :: if_scalar = .false.
     !> Whether to include the buoyancy term in the eddy viscosity.
     logical :: if_buoyancy = .false.
     !> Whether to subtract the horizontal mean of the scalar in the
     !! buoyancy term, as in Abkar & Moin (2017).
     logical :: if_mean_sub = .false.
     !> Name of the scalar field (e.g. temperature).
     character(len=:), allocatable :: scalar_name
     !> The gravity vector, pointing in the direction of gravity.
     real(kind=rp) :: g(3)
     !> Reference temperature, the expansion coefficient is 1 / ref_temp.
     real(kind=rp) :: ref_temp
     !> Cartesian direction (1, 2 or 3) aligned with gravity.
     integer :: gravity_dir = 0
     !> Map from GLL points to levels along gravity, for plane averages.
     type(map_1d_t) :: map_g
     !> SGS diffusivity of the scalar.
     type(field_t), pointer :: alphat => null()
     !> Filter-width tensor (symmetric, 6 components).
     type(field_t) :: g11, g12, g13, g22, g23, g33
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => amd_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          amd_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => amd_free
     !> Compute eddy viscosity (and diffusivity).
     procedure, pass(this) :: compute => amd_compute
  end type amd_t

contains
  !> Constructor.
  !! @param fluid The fluid_scheme_base_t object.
  !! @param json A dictionary with parameters.
  subroutine amd_init(this, fluid, json)
    class(amd_t), intent(inout) :: this
    class(fluid_scheme_base_t), intent(inout), target :: fluid
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: nut_name, alphat_name
    character(len=:), allocatable :: scalar_name
    character(len=:), allocatable :: delta_type
    real(kind=rp), allocatable :: g(:)
    real(kind=rp) :: c, ref_temp
    logical :: if_ext, if_aniso, if_scalar, if_buoyancy, if_mean_sub
    character(len=LOG_SIZE) :: log_buf

    call json_get_or_default(json, "nut_field", nut_name, "nut")
    call json_get_or_default(json, "delta_type", delta_type, "pointwise")
    call json_get_or_default(json, "extrapolation", if_ext, .false.)
    call json_get_or_default(json, "anisotropic", if_aniso, .true.)
    ! Modified Poincare constant, 1/12 for spectrally accurate derivatives.
    call json_get_or_lookup_or_default(json, "c", c, 1.0_rp / 12.0_rp)

    ! Scalar eddy diffusivity (Eq. 18) is computed if a scalar is given.
    if_scalar = json%valid_path("scalar_field")
    if (if_scalar) then
       call json_get(json, "scalar_field", scalar_name)
       call json_get_or_default(json, "alphat_field", alphat_name, &
            "temperature_alphat")
    else
       scalar_name = ""
       alphat_name = ""
    end if

    ! Buoyancy contribution to the eddy viscosity (Eq. 17).
    call json_get_or_default(json, "buoyancy", if_buoyancy, .false.)
    if (if_buoyancy) then
       if (.not. if_scalar) then
          call neko_error("AMD: buoyancy requires scalar_field")
       end if
       call json_get_or_lookup(json, "g", g)
       if (size(g) .ne. 3) then
          call neko_error("AMD: The gravity vector should have 3 components")
       end if
       call json_get_or_lookup(json, "reference_temperature", ref_temp)
       if (ref_temp .le. 0.0_rp) then
          call neko_error("AMD: The reference temperature should be positive")
       end if
       call json_get_or_default(json, "subtract_horizontal_mean", &
            if_mean_sub, .true.)
    else
       g = [0.0_rp, 0.0_rp, 0.0_rp]
       ref_temp = 1.0_rp
       if_mean_sub = .false.
    end if

    call neko_log%section('LES model')
    write(log_buf, '(A)') 'Model : AMD'
    call neko_log%message(log_buf)
    write(log_buf, '(A, A)') 'Delta evaluation : ', delta_type
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'c : ', c
    call neko_log%message(log_buf)
    write(log_buf, '(A, L1)') 'extrapolation : ', if_ext
    call neko_log%message(log_buf)
    write(log_buf, '(A, L1)') 'anisotropic : ', if_aniso
    call neko_log%message(log_buf)
    write(log_buf, '(A, L1)') 'scalar diffusivity : ', if_scalar
    call neko_log%message(log_buf)
    if (if_scalar) then
       write(log_buf, '(A, A)') 'scalar field : ', scalar_name
       call neko_log%message(log_buf)
    end if
    write(log_buf, '(A, L1)') 'buoyancy : ', if_buoyancy
    call neko_log%message(log_buf)
    if (if_buoyancy) then
       write(log_buf, '(A, L1)') 'subtract horizontal mean : ', if_mean_sub
       call neko_log%message(log_buf)
    end if
    call neko_log%end_section()

    call amd_init_from_components(this, fluid, c, nut_name, delta_type, &
         if_ext, if_aniso, if_scalar, scalar_name, alphat_name, &
         if_buoyancy, g, ref_temp, if_mean_sub)

    deallocate(g)

  end subroutine amd_init

  !> Constructor from components.
  !! @param fluid The fluid_scheme_base_t object.
  !! @param c The modified Poincare constant.
  !! @param nut_name The name of the SGS viscosity field.
  !! @param delta_type The type of filter size.
  !! @param if_ext Whether to extrapolate the velocity.
  !! @param if_aniso Whether to use anisotropic filter widths.
  !! @param if_scalar Whether to compute the scalar SGS diffusivity.
  !! @param scalar_name The name of the scalar field.
  !! @param alphat_name The name of the scalar SGS diffusivity field.
  !! @param if_buoyancy Whether to include the buoyancy term.
  !! @param g The gravity vector, pointing in the direction of gravity.
  !! @param ref_temp The reference temperature (beta = 1 / ref_temp).
  !! @param if_mean_sub Whether to subtract the horizontal mean of the scalar.
  subroutine amd_init_from_components(this, fluid, c, nut_name, delta_type, &
       if_ext, if_aniso, if_scalar, scalar_name, alphat_name, if_buoyancy, &
       g, ref_temp, if_mean_sub)
    class(amd_t), intent(inout) :: this
    class(fluid_scheme_base_t), intent(inout), target :: fluid
    real(kind=rp), intent(in) :: c
    character(len=*), intent(in) :: nut_name, delta_type
    logical, intent(in) :: if_ext, if_aniso, if_scalar, if_buoyancy
    logical, intent(in) :: if_mean_sub
    character(len=*), intent(in) :: scalar_name, alphat_name
    real(kind=rp), intent(in) :: g(3)
    real(kind=rp), intent(in) :: ref_temp
    real(kind=rp) :: gmag

    call this%free()
    call this%init_base(fluid, nut_name, delta_type, if_ext)

    this%c = c
    this%if_aniso = if_aniso
    this%if_scalar = if_scalar
    this%if_buoyancy = if_buoyancy
    this%if_mean_sub = if_mean_sub
    this%scalar_name = scalar_name
    this%g = g
    this%ref_temp = ref_temp

    ! The scalar SGS diffusivity is a registry field, so that the scalar
    ! scheme can find it by name (scalar "alphat" -> "alphat_field").
    if (if_scalar) then
       call neko_registry%add_field(fluid%dm_Xh, trim(alphat_name), .true.)
       this%alphat => neko_registry%get_field(trim(alphat_name))
    end if

    ! Filter-width tensor, as in the anisotropic Vreman model.
    call this%g11%init(this%coef%dof, "amd_g11")
    call this%g12%init(this%coef%dof, "amd_g12")
    call this%g13%init(this%coef%dof, "amd_g13")
    call this%g22%init(this%coef%dof, "amd_g22")
    call this%g23%init(this%coef%dof, "amd_g23")
    call this%g33%init(this%coef%dof, "amd_g33")
    if (this%if_aniso) then
       call this%compute_delta_tensor(this%g11, this%g12, this%g13, &
            this%g22, this%g23, this%g33)
    else
       ! Isotropic filter width: G = les_delta**2 * I
       this%g11%x = this%delta%x**2
       this%g22%x = this%delta%x**2
       this%g33%x = this%delta%x**2
       this%g12%x = 0.0_rp
       this%g13%x = 0.0_rp
       this%g23%x = 0.0_rp
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(this%g11%x, this%g11%x_d, this%g11%size(), &
               HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(this%g12%x, this%g12%x_d, this%g12%size(), &
               HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(this%g13%x, this%g13%x_d, this%g13%size(), &
               HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(this%g22%x, this%g22%x_d, this%g22%size(), &
               HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(this%g23%x, this%g23%x_d, this%g23%size(), &
               HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(this%g33%x, this%g33%x_d, this%g33%size(), &
               HOST_TO_DEVICE, sync = .true.)
       end if
    end if

    ! The horizontal mean of the scalar depends on the coordinate along
    ! gravity only, so gravity must be aligned with a Cartesian axis.
    if (if_buoyancy .and. if_mean_sub) then
       gmag = sqrt(sum(g**2))
       if (gmag .le. 0.0_rp) then
          call neko_error("AMD: The gravity vector must be nonzero")
       end if
       this%gravity_dir = maxloc(abs(g), 1)
       if (abs(abs(g(this%gravity_dir)) - gmag) .gt. 1e-6_rp * gmag) then
          call neko_error("AMD: subtract_horizontal_mean requires gravity &
               &aligned with x, y or z")
       end if
       call this%map_g%init(this%coef, this%gravity_dir, 1e-7_rp)
    end if

  end subroutine amd_init_from_components

  !> Destructor for the amd_t class.
  subroutine amd_free(this)
    class(amd_t), intent(inout) :: this

    nullify(this%alphat)
    call this%g11%free()
    call this%g12%free()
    call this%g13%free()
    call this%g22%free()
    call this%g23%free()
    call this%g33%free()
    if (allocated(this%scalar_name)) deallocate(this%scalar_name)
    call this%map_g%free()
    this%gravity_dir = 0

    call this%free_base()
  end subroutine amd_free

  !> Compute eddy viscosity (and diffusivity).
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine amd_compute(this, t, tstep)
    class(amd_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    type(field_t), pointer :: u, v, w, u_e, v_e, w_e

    if (this%if_ext) then
       ! Extrapolate the velocity fields (same as in the other LES models)
       associate(ulag => this%ulag, vlag => this%vlag, &
            wlag => this%wlag, ext_bdf => this%ext_bdf)

         u => neko_registry%get_field_by_name("u")
         v => neko_registry%get_field_by_name("v")
         w => neko_registry%get_field_by_name("w")
         u_e => neko_registry%get_field_by_name("u_e")
         v_e => neko_registry%get_field_by_name("v_e")
         w_e => neko_registry%get_field_by_name("w_e")

         call this%sumab%compute_fluid(u_e, v_e, w_e, u, v, w, &
              ulag, vlag, wlag, ext_bdf%advection_coeffs%x, ext_bdf%nadv)

       end associate
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call amd_compute_device(this%if_ext, t, tstep, this%coef, &
            this%nut, this%alphat, this%c, this%g11, this%g12, this%g13, &
            this%g22, this%g23, this%g33, this%if_scalar, this%scalar_name, &
            this%if_buoyancy, this%g, this%ref_temp, this%if_mean_sub, &
            this%gravity_dir, this%map_g)
    else
       call amd_compute_cpu(this%if_ext, t, tstep, this%coef, &
            this%nut, this%alphat, this%c, this%g11, this%g12, this%g13, &
            this%g22, this%g23, this%g33, this%if_scalar, this%scalar_name, &
            this%if_buoyancy, this%g, this%ref_temp, this%if_mean_sub, &
            this%gravity_dir, this%map_g)
    end if

  end subroutine amd_compute

end module amd
