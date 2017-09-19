module ModIE_Interface

  use ModKind

  real(real8_), allocatable, dimension(:,:,:) :: IEr3_HaveLats, IEr3_HaveMLTs
  real(real8_), allocatable, dimension(:,:,:) :: IEr3_HavePotential
  real(real8_), allocatable, dimension(:,:,:) :: IEr3_HaveEFlux
  real(real8_), allocatable, dimension(:,:,:) :: IEr3_HaveAveE

  real (kind=real8_) :: IEd_CurrentTime
  integer             :: IEi_HavenLats
  integer             :: IEi_HavenMLTs
  integer             :: IEi_HavenBLKs
  integer             :: IEi_HavenTimes

  real (kind=real8_)               :: GMd_NeedTime = -1.0e32
  real(real8_), allocatable, dimension(:,:) :: GMr2_NeedLats, GMr2_NeedMLTs
  real(real8_), allocatable, dimension(:,:) :: GMr2_NeedPotential
  real(real8_), allocatable, dimension(:,:) :: GMr2_NeedEFlux
  real(real8_), allocatable, dimension(:,:) :: GMr2_NeedAveE
  integer                           :: GMi_NeednLats
  integer                           :: GMi_NeednMLTs
  integer                           :: GMi_NeednTimes
  integer, allocatable, dimension(:,:,:) :: GMi3_InterpolationIndices
  real(real8_), allocatable, dimension(:,:,:)    :: GMr3_InterpolationRatios

  real (kind=real8_)               :: IMd_NeedTime = -1.0e32
  real(real8_), allocatable, dimension(:,:) :: IMr2_NeedLats, IMr2_NeedMLTs
  real(real8_), allocatable, dimension(:,:) :: IMr2_NeedPotential
  real(real8_), allocatable, dimension(:,:) :: IMr2_NeedEFlux
  real(real8_), allocatable, dimension(:,:) :: IMr2_NeedAveE
  integer                           :: IMi_NeednLats
  integer                           :: IMi_NeednMLTs
  integer                           :: IMi_NeednTimes
  integer, allocatable, dimension(:,:,:) :: IMi3_InterpolationIndices
  real(real8_), allocatable, dimension(:,:,:)    :: IMr3_InterpolationRatios

  real (kind=real8_)               :: IOd_NeedTime = -1.0e32
  real(real8_), allocatable, dimension(:,:) :: IOr2_NeedLats, IOr2_NeedMLTs
  real(real8_), allocatable, dimension(:,:) :: IOr2_NeedPotential
  real(real8_), allocatable, dimension(:,:) :: IOr2_NeedEFlux
  real(real8_), allocatable, dimension(:,:) :: IOr2_NeedAveE
  integer                           :: IOi_NeednLats
  integer                           :: IOi_NeednMLTs
  integer                           :: IOi_NeednTimes
  integer, allocatable, dimension(:,:,:) :: IOi3_InterpolationIndices
  real(real8_), allocatable, dimension(:,:,:)    :: IOr3_InterpolationRatios
  real(real8_) :: IOr_NeedIMFBz   = -1.0e32
  real(real8_) :: IOr_NeedIMFBy   = -1.0e32 
  real(real8_) :: IOr_NeedSWV     = -1.0e32 
  real(real8_) :: IOr_NeedHPI     = -1.0e32 
  real(real8_) :: IOr_NeedHPINorm = -1.0e32 
  real(real8_) :: IOr_NeedKp      = -1.0e32 
  logical :: IOl_IsNorth  = .true.

  integer, parameter                :: IE_Closest_     = 1
  integer, parameter                :: IE_After_       = 2
  integer, parameter                :: IE_Interpolate_ = 3

  character (len=100) :: IE_NameOfEFieldModel
  character (len=100) :: IE_NameOfAuroralModel
  character (len=100) :: IE_NameOfSolarModel
  character (len=100) :: IE_NameOfModelDir

  logical :: UseGridBasedIE
  logical :: UseAMIE
  logical :: UseSPS
  logical :: UsePPS

end module ModIE_Interface
