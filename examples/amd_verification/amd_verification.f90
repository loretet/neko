! Verification of the AMD model on frozen analytic fields (see README.md).
! u_i = A_ij x_j (constant velocity gradient A),
! theta = 300 + 6 x - 4 y + 0.05 z + 0.3 z^2.
! The SEM derivatives of these polynomials are exact, so nut and alphat
! computed by Neko must match the analytic AMD formulas to round-off.
module user
  use neko
  implicit none

  real(kind=rp) :: A(3,3)
  logical :: mean_sub = .true.
  character(len=16) :: mode
  real(kind=rp), parameter :: CPOIN = 1.0_rp / 12.0_rp
  real(kind=rp) :: BG(3) = [0.0_rp, 0.0_rp, 9.81_rp / 300.0_rp]

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    ! A(i,j) = du_i/dx_j, selected by the AMD_TEST_MODE environment variable
    call get_environment_variable("AMD_TEST_MODE", mode)
    A = reshape([0.1_rp, 0.3_rp, -0.2_rp, 0.5_rp, -0.4_rp, 0.7_rp, &
         0.2_rp, 0.6_rp, 0.3_rp], [3, 3])
    select case (trim(mode))
    case ("shear")      ! u = 0.8 z: laminar shear, AMD must give nut = 0
       A = 0.0_rp; A(1,3) = 0.8_rp
    case ("rotation")   ! solid-body rotation about z: S = 0
       A = 0.0_rp; A(1,2) = -0.5_rp; A(2,1) = 0.5_rp
    case ("negA")       ! reversed flow: production changes sign -> clipping
       A = -A
    case ("nomean")
       mean_sub = .false.
    case ("neutral")    ! no buoyancy term (run amd_verification_neutral.case)
       BG = 0.0_rp
    end select
    if (trim(mode) .eq. "skew") user%mesh_setup => skew_mesh
    user%initial_conditions => user_ic
    user%compute => user_check
  end subroutine user_setup

  ! Shear the mesh, x -> x + 0.3 y + 0.2 z, so that the elements are no
  ! longer aligned with the axes and the filter-width tensor G has
  ! off-diagonal components. Planes of constant z are unchanged.
  subroutine skew_mesh(msh, time)
    type(mesh_t), intent(inout) :: msh
    type(time_state_t), intent(in) :: time
    integer :: i

    do i = 1, size(msh%points)
       msh%points(i)%x(1) = msh%points(i)%x(1) + 0.3_rp * msh%points(i)%x(2) &
            + 0.2_rp * msh%points(i)%x(3)
    end do
  end subroutine skew_mesh

  subroutine user_ic(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: u, v, w, p, s
    real(kind=rp) :: x(3)
    integer :: i

    if (scheme_name .eq. 'fluid') then
       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")
       p => fields%get("p")
       call field_rzero(p)
       do i = 1, u%dof%size()
          x = [u%dof%x%x(i,1,1,1), u%dof%y%x(i,1,1,1), u%dof%z%x(i,1,1,1)]
          u%x(i,1,1,1) = sum(A(1,:) * x)
          v%x(i,1,1,1) = sum(A(2,:) * x)
          w%x(i,1,1,1) = sum(A(3,:) * x)
       end do
    else
       s => fields%get(scheme_name)
       do i = 1, s%dof%size()
          x = [s%dof%x%x(i,1,1,1), s%dof%y%x(i,1,1,1), s%dof%z%x(i,1,1,1)]
          s%x(i,1,1,1) = 300.0_rp + 6.0_rp * x(1) - 4.0_rp * x(2) &
               + 0.05_rp * x(3) + 0.3_rp * x(3)**2
       end do
    end if
  end subroutine user_ic

  ! Analytic reference written with matrix operations, independently of the
  ! loop structure used in amd_cpu.f90.
  subroutine user_check(time)
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: nut, alphat
    real(kind=rp) :: S(3,3), D(3,3), G(3,3), gt(3), gtb(3), bu(3), q(3)
    real(kind=rp) :: nut_ref, alphat_ref, err_nut, err_alphat
    real(kind=rp) :: max_nut, max_alphat
    real(kind=rp) :: z
    integer :: i, nclip, lx, ii, jj, kk, e
    character(len=LOG_SIZE) :: log_buf

    if (time%tstep .lt. 2) return
    if (time%tstep .gt. 2) return

    nut => neko_registry%get_field("nut")
    alphat => neko_registry%get_field("temperature_alphat")
    lx = nut%Xh%lx

    S = 0.5_rp * (A + transpose(A))
    err_nut = 0.0_rp; err_alphat = 0.0_rp
    max_nut = 0.0_rp; max_alphat = 0.0_rp
    nclip = 0
    do e = 1, nut%msh%nelv
    do kk = 1, lx
    do jj = 1, lx
    do ii = 1, lx
       i = ii + lx * (jj - 1) + lx * lx * (kk - 1) + lx**3 * (e - 1)
       z = nut%dof%z%x(i,1,1,1)
       ! Filter-width tensor from the GLL grid steps ("pointwise")
       G = CPOIN * filter_width_tensor(nut%dof, ii, jj, kk, e, lx)
       ! D_ij = (dk^ u_i)(dk^ u_j) = A G A^T
       D = matmul(A, matmul(G, transpose(A)))
       gt = [6.0_rp, -4.0_rp, 0.05_rp + 0.6_rp * z]
       gtb = [6.0_rp, -4.0_rp, 0.0_rp]      ! gradient of theta - <theta>
       if (.not. mean_sub) gtb = gt
       bu = matmul(BG, A)                    ! bu_k = BG_i A_ik
       nut_ref = max(0.0_rp, (-sum(D * S) + dot_product(bu, matmul(G, gtb))) &
            / sum(A**2))
       q = matmul(gt, A)                     ! q_k = A_ik gt_i
       alphat_ref = max(0.0_rp, -dot_product(q, matmul(G, gt)) / sum(gt**2))
       if (alphat_ref .eq. 0.0_rp) nclip = nclip + 1

       err_nut = max(err_nut, abs(nut%x(i,1,1,1) - nut_ref))
       err_alphat = max(err_alphat, abs(alphat%x(i,1,1,1) - alphat_ref))
       max_nut = max(max_nut, nut_ref)
       max_alphat = max(max_alphat, alphat_ref)
    end do
    end do
    end do
    end do
    err_nut = glmax([err_nut], 1)
    err_alphat = glmax([err_alphat], 1)
    max_nut = glmax([max_nut], 1)
    max_alphat = glmax([max_alphat], 1)
    nclip = glimax([nclip], 1)

    write(log_buf, '(A,ES12.4,A,ES12.4)') 'AMD TEST nut    max ', max_nut, &
         '  max abs err ', err_nut
    call neko_log%message(log_buf)
    write(log_buf, '(A,ES12.4,A,ES12.4)') 'AMD TEST alphat max ', max_alphat, &
         '  max abs err ', err_alphat
    call neko_log%message(log_buf)
    write(log_buf, '(A,A,A,I8)') 'AMD TEST mode ', trim(mode), &
         ' clipped alphat points (max rank) ', nclip
    call neko_log%message(log_buf)
  end subroutine user_check

  ! G = dr dr^T + ds ds^T + dt dt^T, with dr, ds, dt the grid steps around
  ! GLL point (i,j,k) of element e along the local directions r, s, t.
  function filter_width_tensor(dof, i, j, k, e, lx) result(G)
    type(dofmap_t), intent(in) :: dof
    integer, intent(in) :: i, j, k, e, lx
    real(kind=rp) :: G(3,3), dr(3), ds(3), dt(3)
    integer :: m

    dr = point(min(i + 1, lx), j, k) - point(max(i - 1, 1), j, k)
    dr = dr / real(min(i + 1, lx) - max(i - 1, 1), rp)
    ds = point(i, min(j + 1, lx), k) - point(i, max(j - 1, 1), k)
    ds = ds / real(min(j + 1, lx) - max(j - 1, 1), rp)
    dt = point(i, j, min(k + 1, lx)) - point(i, j, max(k - 1, 1))
    dt = dt / real(min(k + 1, lx) - max(k - 1, 1), rp)
    do m = 1, 3
       G(:,m) = dr * dr(m) + ds * ds(m) + dt * dt(m)
    end do
  contains
    function point(ip, jp, kp) result(x)
      integer, intent(in) :: ip, jp, kp
      real(kind=rp) :: x(3)
      x = [dof%x%x(ip, jp, kp, e), dof%y%x(ip, jp, kp, e), &
           dof%z%x(ip, jp, kp, e)]
    end function point
  end function filter_width_tensor

end module user
