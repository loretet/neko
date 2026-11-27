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
!> Implements the device kernel for the `amd_t` type.
!! @note Not implemented yet: the device version should follow
!! `vreman_device.f90` (scratch fields on the device and a pointwise kernel
!! for each backend: CUDA, HIP, OpenCL and Metal).
module amd_device
  use num_types, only : rp
  use field, only : field_t
  use coefs, only : coef_t
  use map_1d, only : map_1d_t
  use utils, only : neko_error
  implicit none
  private

  public :: amd_compute_device

contains

  !> Compute the AMD eddy viscosity on the device (not implemented yet).
  !! The arguments are the same as for `amd_compute_cpu`.
  subroutine amd_compute_device(if_ext, t, tstep, coef, nut, alphat, c, &
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

    call neko_error("AMD: the device (GPU) backend is not implemented yet")

  end subroutine amd_compute_device

end module amd_device
