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
!> Implements the CPU kernel for the `amd_t` type.
!!
!! =====================================================================
!!  STUDENT VERSION.  Fill in every block marked  TODO.
!! =====================================================================
!!  - The file compiles and runs as it is: every unfinished quantity is
!!    set to zero, so you can build and test after each step.
!!  - Work through the TODOs in order (TODO 1, 2, ...). After each stage,
!!    rebuild Neko, rerun makeneko and run the verification case (see the
!!    guide): the reported errors must go down to round-off (~1e-15).
!!  - Notation used in the comments:
!!      a_ij   = d u_i / d x_j       (velocity gradient tensor, matrix A)
!!      S_ij   = (a_ij + a_ji) / 2   (strain-rate tensor)
!!      G_ab   = filter-width tensor (computed in les_model.f90); for an
!!               axis-aligned element G = diag(dx**2, dy**2, dz**2)
!!      (dk^ f)(dk^ g) = c * (df/dx_a) G_ab (dg/dx_b)
!!               (product of two scaled gradients, see the guide)
!!      s      = the scalar (temperature), s' = s - <s> (horizontal mean)
!!    Repeated indices are summed (i, j, k, a, b = 1, 2, 3).
!! =====================================================================
module amd_cpu
  use num_types, only : rp
  use math, only : NEKO_EPS, col2
  use scratch_registry, only : neko_scratch_registry
  use registry, only : neko_registry
  use field, only : field_t
  use field_list, only : field_list_t
  use matrix, only : matrix_t
  use map_1d, only : map_1d_t
  use operators, only : dudxyz, grad
  use coefs, only : coef_t
  use gs_ops, only : GS_OP_ADD
  implicit none
  private

  public :: amd_compute_cpu

contains

  !> Compute the AMD eddy viscosity (and scalar eddy diffusivity) on the CPU.
  !! Implements Eqs. (17)-(18) of Abkar & Moin (2017):
  !!   nut    = max(0, [-(dk^ u_i)(dk^ u_j) S_ij
  !!                    - g_i / T_ref (dk^ u_i)(dk^ s')] / (dl u_m dl u_m))
  !!   alphat = max(0, -(dk^ u_i)(dk^ s) di s / (dl s dl s))
  !! The products of scaled gradients are evaluated with the filter-width
  !! tensor G: (dk^ f)(dk^ g) = c * (df/dx_a) G_ab (dg/dx_b).
  !! For an axis-aligned element G = diag(dx^2, dy^2, dz^2), which gives
  !! the usual sum_k c delta_k^2 (df/dx_k)(dg/dx_k).
  !! @param if_ext If extrapolate the velocity field to evaluate
  !! @param t The time value.
  !! @param tstep The current time-step.
  !! @param coef SEM coefficients.
  !! @param nut The SGS viscosity field (output).
  !! @param alphat The SGS scalar diffusivity field (output, if if_scalar).
  !! @param c The modified Poincare constant.
  !! @param g11, g12, g13, g22, g23, g33 The filter-width tensor.
  !! @param if_scalar Whether to compute alphat.
  !! @param scalar_name The name of the scalar field.
  !! @param if_buoyancy Whether to include the buoyancy term in nut.
  !! @param g The gravity vector.
  !! @param ref_temp The reference temperature.
  !! @param if_mean_sub Whether to subtract the horizontal mean of the scalar.
  !! @param gravity_dir The Cartesian direction aligned with gravity.
  !! @param map_g Map from GLL points to levels along gravity.
  subroutine amd_compute_cpu(if_ext, t, tstep, coef, nut, alphat, c, &
       g11, g12, g13, g22, g23, g33, if_scalar, scalar_name, if_buoyancy, &
       g, ref_temp, if_mean_sub, gravity_dir, map_g)
    logical, intent(in) :: if_ext
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(in) :: coef
    type(field_t), intent(inout) :: nut
    type(field_t), pointer, intent(inout) :: alphat
    real(kind=rp), intent(in) :: c
    type(field_t), intent(in) :: g11, g12, g13, g22, g23, g33
    logical, intent(in) :: if_scalar, if_buoyancy, if_mean_sub
    character(len=*), intent(in) :: scalar_name
    real(kind=rp), intent(in) :: g(3)
    real(kind=rp), intent(in) :: ref_temp
    integer, intent(in) :: gravity_dir
    type(map_1d_t), intent(inout) :: map_g

    ! Velocity gradient tensor, a_ij = d u_i / d x_j
    type(field_t), pointer :: a11, a12, a13, a21, a22, a23, a31, a32, a33
    type(field_t), pointer :: u, v, w
    ! Scalar and its gradient
    type(field_t), pointer :: s, dsdx, dsdy, dsdz
    integer :: temp_indices(9), temp_indices_s(3)
    integer :: e, i, j, k
    ! Pointwise work variables
    real(kind=rp) :: a(3,3), gt(3,3), h(3,3), d(3,3), sr(3,3)
    real(kind=rp) :: shear_prod, buoy_prod, aijaij
    real(kind=rp) :: gs(3), gs_b(3), bu(3), q(3), grad_s_sq
    ! Plane averages of the scalar gradient along gravity
    real(kind=rp), allocatable :: mean_grad(:)

    if (if_ext) then
       u => neko_registry%get_field_by_name("u_e")
       v => neko_registry%get_field_by_name("v_e")
       w => neko_registry%get_field_by_name("w_e")
    else
       u => neko_registry%get_field_by_name("u")
       v => neko_registry%get_field_by_name("v")
       w => neko_registry%get_field_by_name("w")
    end if

    call neko_scratch_registry%request_field(a11, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(a12, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(a13, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(a21, temp_indices(4), .false.)
    call neko_scratch_registry%request_field(a22, temp_indices(5), .false.)
    call neko_scratch_registry%request_field(a23, temp_indices(6), .false.)
    call neko_scratch_registry%request_field(a31, temp_indices(7), .false.)
    call neko_scratch_registry%request_field(a32, temp_indices(8), .false.)
    call neko_scratch_registry%request_field(a33, temp_indices(9), .false.)

    ! ------------------------------------------------------------------
    ! TODO 1: velocity gradient tensor a_ij = d u_i / d x_j.
    !
    ! dudxyz(result, f, dr/dx_j, ds/dx_j, dt/dx_j, coef) computes d f / d x_j
    ! (chain rule from the element coordinates r, s, t to x_j).
    ! The first component is given. Write the other eight:
    !   - the field f is u (row 1), v (row 2) or w (row 3);
    !   - the metric terms select the direction x_j:
    !       x: coef%drdx, coef%dsdx, coef%dtdx
    !       y: coef%drdy, coef%dsdy, coef%dtdy
    !       z: coef%drdz, coef%dsdz, coef%dtdz
    ! Compare with vreman_cpu.f90, where the same tensor is computed.
    ! ------------------------------------------------------------------
    call dudxyz (a11%x, u%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    ! call dudxyz (a12%x, ...)        ! TODO: du/dy
    ! call dudxyz (a13%x, ...)        ! TODO: du/dz
    ! call dudxyz (a21%x, ...)        ! TODO: dv/dx
    ! call dudxyz (a22%x, ...)        ! TODO: dv/dy
    ! call dudxyz (a23%x, ...)        ! TODO: dv/dz
    ! call dudxyz (a31%x, ...)        ! TODO: dw/dx
    ! call dudxyz (a32%x, ...)        ! TODO: dw/dy
    ! call dudxyz (a33%x, ...)        ! TODO: dw/dz
    ! (remove the next line when all eight calls above are written)
    call zero_fields_tmp(a12, a13, a21, a22, a23, a31, a32, a33)

    ! ------------------------------------------------------------------
    ! TODO 2: make the gradients continuous across elements.
    !
    ! dudxyz works element by element, so at a GLL point shared by several
    ! elements each element gives a slightly different derivative.
    ! average_across_elements (bottom of this file) sums the values of all
    ! elements sharing the point and divides by their number.
    ! (Vreman only sums them, which works because its formula does not
    ! change when all gradients are scaled; AMD mixes velocity and scalar
    ! gradients, so a proper average is needed.)
    ! The first one is given; do the same for the other eight components.
    ! ------------------------------------------------------------------
    call average_across_elements(a11, coef)
    ! TODO: a12, a13, a21, a22, a23, a31, a32, a33

    ! Scalar gradient (given). grad() calls dudxyz three times.
    if (if_scalar) then
       s => neko_registry%get_field_by_name(scalar_name)
       call neko_scratch_registry%request_field(dsdx, temp_indices_s(1), &
            .false.)
       call neko_scratch_registry%request_field(dsdy, temp_indices_s(2), &
            .false.)
       call neko_scratch_registry%request_field(dsdz, temp_indices_s(3), &
            .false.)
       call grad(dsdx%x, dsdy%x, dsdz%x, s%x, coef)
       call average_across_elements(dsdx, coef)
       call average_across_elements(dsdy, coef)
       call average_across_elements(dsdz, coef)
    end if

    ! Horizontal mean of the scalar gradient along gravity, d<s>/dz (given).
    ! mean_grad(l) holds <d s/dz> on GLL level l; map_g%pt_lvl(i,1,1,e)
    ! gives the level of point i in element e.
    ! Since <s> depends only on the coordinate along gravity, the gradient
    ! of the fluctuation s' = s - <s> differs from grad(s) only in that
    ! component: d s'/dz = d s/dz - <d s/dz>.
    allocate(mean_grad(max(1, map_g%n_gll_lvls)))
    mean_grad = 0.0_rp
    if (if_buoyancy .and. if_mean_sub) then
       select case (gravity_dir)
       case (1)
          call plane_average(mean_grad, dsdx, map_g)
       case (2)
          call plane_average(mean_grad, dsdy, map_g)
       case (3)
          call plane_average(mean_grad, dsdz, map_g)
       end select
    end if

    gs = 0.0_rp
    gs_b = 0.0_rp

    !$omp parallel do private(e, i, j, k, a, gt, h, d, sr, shear_prod, &
    !$omp& buoy_prod, aijaij, gs, gs_b, bu, q, grad_s_sq)
    do e = 1, coef%msh%nelv
       do i = 1, coef%Xh%lxyz
          ! ------------------------------------------------------------
          ! TODO 3: copy the velocity gradient at this point into the 3x3
          ! matrix a, with a(i,j) = d u_i / d x_j. Row 1 is given.
          ! ------------------------------------------------------------
          a = 0.0_rp
          a(1,1) = a11%x(i,1,1,e)
          a(1,2) = a12%x(i,1,1,e)
          a(1,3) = a13%x(i,1,1,e)
          ! TODO: rows 2 (v) and 3 (w)

          ! ------------------------------------------------------------
          ! TODO 4: the scaled filter-width tensor gt = c * G.
          ! G is symmetric and only 6 components are stored (g11, g12,
          ! g13, g22, g23, g33); fill the upper triangle from them and copy
          ! it to the lower triangle. gt(1,1) is given.
          ! ------------------------------------------------------------
          gt = 0.0_rp
          gt(1,1) = c * g11%x(i,1,1,e)
          ! TODO: gt(1,2), gt(1,3), gt(2,2), gt(2,3), gt(3,3)
          ! TODO: gt(2,1), gt(3,1), gt(3,2) (symmetry)

          ! ------------------------------------------------------------
          ! TODO 5: strain-rate tensor S = (A + A^T) / 2 (one line, using
          ! the Fortran intrinsic transpose).
          ! ------------------------------------------------------------
          sr = 0.0_rp   ! TODO

          ! ------------------------------------------------------------
          ! TODO 6: the tensor D_ij = (dk^ u_i)(dk^ u_j), i.e. D = A (c G) A^T.
          ! It is computed in two steps, as in the anisotropic Vreman model:
          !   h(k,j) = sum_b gt(k,b) * a(j,b)     (h = (c G) A^T)
          !   d(j,k) = sum_m a(j,m) * h(m,k)      (D = A h)
          ! The loops are given: write the two sums (3 terms each).
          ! ------------------------------------------------------------
          do j = 1, 3
             do k = 1, 3
                h(k,j) = 0.0_rp   ! TODO
             end do
          end do

          do j = 1, 3
             do k = 1, 3
                d(j,k) = 0.0_rp   ! TODO
             end do
          end do

          ! ------------------------------------------------------------
          ! TODO 7: shear production  -D_ij S_ij  (double contraction: the
          ! sum of the element-wise products of the two matrices) and the
          ! denominator a_lm a_lm (sum of the squares of all nine a_ij).
          ! Hint: the intrinsic sum() applied to a whole array.
          ! ------------------------------------------------------------
          shear_prod = 0.0_rp   ! TODO
          aijaij = 0.0_rp       ! TODO

          buoy_prod = 0.0_rp
          if (if_scalar) then
             ! Scalar gradient
             gs(1) = dsdx%x(i,1,1,e)
             gs(2) = dsdy%x(i,1,1,e)
             gs(3) = dsdz%x(i,1,1,e)

             if (if_buoyancy) then
                ! ------------------------------------------------------
                ! TODO 8: bu_k = -g_i * a_ik / ref_temp  (sum over i).
                ! This is d/dx_k of the velocity component along gravity,
                ! times the buoyancy parameter 1 / ref_temp. Do not forget
                ! the minus: g points in the direction of gravity (it is
                ! negative in the input), so for g = (0, 0, -9.81) you
                ! should get +(9.81 / ref_temp) * d w / d x_k.
                ! Careful: the sum runs over the FIRST index of a.
                ! ------------------------------------------------------
                do k = 1, 3
                   bu(k) = 0.0_rp   ! TODO
                end do

                ! ------------------------------------------------------
                ! TODO 9: gradient of the fluctuation s' = s - <s>.
                ! Start from gs and, if if_mean_sub is true, subtract the
                ! plane-averaged gradient from the component along gravity:
                !   gs_b(gravity_dir) = gs_b(gravity_dir)
                !                       - mean_grad(map_g%pt_lvl(i,1,1,e))
                ! ------------------------------------------------------
                gs_b = gs
                ! TODO: mean subtraction (inside an  if (if_mean_sub)  block)

                ! ------------------------------------------------------
                ! TODO 10: buoyancy production
                !   -g_i / T_ref (dk^ u_i)(dk^ s') = bu_a gt_ab gs_b(b)
                ! Hint: matmul (matrix times vector) and dot_product.
                ! ------------------------------------------------------
                buoy_prod = 0.0_rp   ! TODO
             end if

             ! ---------------------------------------------------------
             ! TODO 11: q_k = a_ik * (d s / d x_i)  (sum over i): how the
             ! velocity gradient stretches the scalar gradient. Also
             ! compute |grad s|^2.
             ! ---------------------------------------------------------
             do k = 1, 3
                q(k) = 0.0_rp   ! TODO
             end do
             grad_s_sq = 0.0_rp   ! TODO

             ! ---------------------------------------------------------
             ! TODO 12: eddy diffusivity, Eq. (18):
             !   alphat = max(0, - q_a gt_ab gs_b / (|grad s|^2 + eps))
             ! Use NEKO_EPS as eps, to avoid dividing by zero where the
             ! scalar is uniform.
             ! ---------------------------------------------------------
             alphat%x(i,1,1,e) = 0.0_rp   ! TODO
          end if

          ! ------------------------------------------------------------
          ! TODO 13: eddy viscosity, Eq. (17):
          !   nut = max(0, (shear production + buoyancy production)
          !                / (a_lm a_lm + eps))
          ! The max(0, .) is the clipping: a negative value would be
          ! backscatter, which is not allowed in the model.
          ! ------------------------------------------------------------
          nut%x(i,1,1,e) = 0.0_rp   ! TODO
       end do
    end do
    !$omp end parallel do

    ! The filter widths can differ between elements at shared GLL points,
    ! so make the outputs continuous, as in the other LES models.
    call average_across_elements(nut, coef)
    if (if_scalar) then
       call average_across_elements(alphat, coef)
       call neko_scratch_registry%relinquish_field(temp_indices_s)
    end if

    deallocate(mean_grad)
    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine amd_compute_cpu

  !> Average a field across elements: gather-scatter sum, then divide by
  !! the multiplicity of each GLL point.
  !! @param f The field to average.
  !! @param coef SEM coefficients.
  subroutine average_across_elements(f, coef)
    type(field_t), intent(inout) :: f
    type(coef_t), intent(in) :: coef

    call coef%gs_h%op(f, GS_OP_ADD)
    call col2(f%x, coef%mult, f%dof%size())
  end subroutine average_across_elements

  !> Average a field over planes normal to the direction of `map`.
  !! @param mean Output: plane average for each GLL level.
  !! @param f The field to average.
  !! @param map The map from GLL points to levels.
  subroutine plane_average(mean, f, map)
    real(kind=rp), intent(inout) :: mean(:)
    type(field_t), target, intent(in) :: f
    type(map_1d_t), intent(inout) :: map
    type(field_list_t) :: fl
    type(matrix_t) :: avg

    call fl%init(1)
    call fl%assign_to_field(1, f)
    call map%average_planes_fld_lst(avg, fl)
    ! Column 1 holds the coordinates, column 2 the averages.
    mean(1:map%n_gll_lvls) = avg%x(1:map%n_gll_lvls, 2)
    call avg%free()
    call fl%free()
  end subroutine plane_average

  !> Temporary helper for the student version: sets the gradient fields that
  !! are not computed yet to zero. Delete it (and its call) after TODO 1.
  subroutine zero_fields_tmp(f1, f2, f3, f4, f5, f6, f7, f8)
    type(field_t), intent(inout) :: f1, f2, f3, f4, f5, f6, f7, f8
    f1%x = 0.0_rp; f2%x = 0.0_rp; f3%x = 0.0_rp; f4%x = 0.0_rp
    f5%x = 0.0_rp; f6%x = 0.0_rp; f7%x = 0.0_rp; f8%x = 0.0_rp
  end subroutine zero_fields_tmp

end module amd_cpu