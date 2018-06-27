!--------------------------------------------------
!>Specify the KIND value
!--------------------------------------------------
module NumberKinds
    implicit none

    integer, parameter                                  :: KREAL = kind(0.d0)
    integer, parameter                                  :: KINT = kind(1)
end module NumberKinds

!--------------------------------------------------
!>Define the global constant variables
!--------------------------------------------------
module ConstantVariables
    use NumberKinds
    implicit none

    real(KREAL), parameter                              :: PI = 4.0*atan(1.0) !Pi
    real(KREAL), parameter                              :: SMV = tiny(0.0) !Small value to avoid 0/0
    real(KREAL), parameter                              :: UP = 1.0 !Used in sign() function
end module ConstantVariables

!--------------------------------------------------
!>Define the mesh data structure
!--------------------------------------------------
module Mesh
    use NumberKinds
    implicit none

    !--------------------------------------------------
    !Basic derived type
    !--------------------------------------------------
    !Cell center
    type CellCenter
        !Geometry
        real(KREAL)                                     :: x !Cell center coordinates
        real(KREAL)                                     :: length !Cell length
        !Flow field
        real(KREAL)                                     :: conVars(3) !Conservative variables at cell center: density, x-momentum, total energy
        real(KREAL), allocatable, dimension(:)          :: h,b !Distribution function
        real(KREAL), allocatable, dimension(:)          :: sh,sb !Slope of distribution function
    end type CellCenter

    !Cell interface
    type CellInterface
        !Flux
        real(KREAL)                                     :: flux(3) !Conservative variables flux at cell interface
        real(KREAL), allocatable, dimension(:)          :: flux_h,flux_b !Flux of distribution function
    end type CellInterface

    !Index method
    ! --------------------
    !  (i) |  (i)  | (i+1)
    ! face | cell  | face
    ! --------------------
end module Mesh

module ControlParameters
    use ConstantVariables
    use Mesh
    implicit none

    !--------------------------------------------------
    !Variables to control the simulation
    !--------------------------------------------------
    real(KREAL), parameter                              :: CFL = 0.95 !CFL number
    real(KREAL)                                         :: simTime = 0.0 !Current simulation time
    real(KREAL), parameter                              :: MAX_TIME = 250.0 !Maximal simulation time
    integer(KINT), parameter                            :: MAX_ITER = 5E5 !Maximal iteration number
    integer(KINT)                                       :: iter = 1 !Number of iteration
    real(KREAL)                                         :: dt !Global time step

    !Output control
    character(len=24), parameter                        :: HSTFILENAME = "StationaryShockStructure" !History file name
    character(len=24), parameter                        :: RSTFILENAME = "StationaryShockStructure" !Result file name
    integer(KINT), parameter                            :: HSTFILE = 20 !History file ID
    integer(KINT), parameter                            :: RSTFILE = 21 !Result file ID

    !Gas propeties
    integer(KINT), parameter                            :: CK = 2 !Internal degree of freedom, here 2 denotes monatomic gas
    real(KREAL), parameter                              :: GAMMA = real(CK+3,KREAL)/real(CK+1,KREAL) !Ratio of specific heat
    real(KREAL), parameter                              :: OMEGA = 0.72 !Temperature dependence index in HS/VHS/VSS model
    real(KREAL), parameter                              :: PR = 2.0/3.0 !Prandtl number
    real(KREAL), parameter                              :: KN = 1.0 !Knudsen number in reference state
    real(KREAL), parameter                              :: ALPHA_REF = 1.0 !Coefficient in HS model
    real(KREAL), parameter                              :: OMEGA_REF = 0.5 !Coefficient in HS model
    real(KREAL), parameter                              :: MU_REF = 5.0*(ALPHA_REF+1.0)*(ALPHA_REF+2.0)*sqrt(PI)/(4.0*ALPHA_REF*(5.0-2.0*OMEGA_REF)*(7.0-2.0*OMEGA_REF))*KN !Viscosity coefficient in reference state
    real(KREAL), parameter                              :: MA = 8.0 !Mach number
    !Geometry
    real(KREAL), parameter                              :: START_POINT = 0.0, END_POINT = 50.0
    integer(KINT), parameter                            :: POINTS_NUM = 100
    integer(KINT), parameter                            :: IXMIN = 1 , IXMAX = POINTS_NUM !Cell index range
    integer(KINT), parameter                            :: GHOST_NUM = 1 !Ghost cell number

    !--------------------------------------------------
    !Initial flow field
    !--------------------------------------------------
    !Index method
    !-------------------------------
    !| (i-1) |  (i) |  (i) | (i+1) |
    !| cell  | face | cell | face  |
    !-------------------------------
    type(CellCenter)                                    :: ctr(IXMIN-GHOST_NUM:IXMAX+GHOST_NUM) !Cell center (with ghost cell)
    type(CellInterface)                                 :: vface(IXMIN-GHOST_NUM+1:IXMAX+GHOST_NUM) !Vertical cell interfaces

    !Initial value
    !Upstream condition (before shock)
    real(KREAL), parameter                              :: PRIM_L1 = 1.0 !Density
    real(KREAL), parameter                              :: PRIM_L2 = MA*sqrt(GAMMA/2.0) !X-velocity
    real(KREAL), parameter                              :: PRIM_L3 = 1.0 !Lambda=1/temperature
    real(KREAL), parameter, dimension(3)                :: PRIM_LEFT = [PRIM_L1,PRIM_L2,PRIM_L3]

    !Downstream condition (after shock)
    real(KREAL), parameter                              :: MA_R = sqrt((MA**2*(GAMMA-1)+2)/(2*GAMMA*MA**2-(GAMMA-1)))
    real(KREAL), parameter                              :: RATIO_T = (1+(GAMMA-1)/2*MA**2)*(2*GAMMA/(GAMMA-1)*MA**2-1)/(MA**2*(2*GAMMA/(GAMMA-1)+(GAMMA-1)/2))

    real(KREAL), parameter                              :: PRIM_R1 = PRIM_LEFT(1)*(GAMMA+1)*MA**2/((GAMMA-1)*MA**2+2)
    real(KREAL), parameter                              :: PRIM_R2 = MA_R*sqrt(GAMMA/2.0)*sqrt(RATIO_T)
    real(KREAL), parameter                              :: PRIM_R3 = PRIM_LEFT(3)/RATIO_T
    real(KREAL), parameter, dimension(3)                :: PRIM_RIGHT = [PRIM_R1,PRIM_R2,PRIM_R3]

    !--------------------------------------------------
    !Discrete velocity space
    !--------------------------------------------------
    integer(KINT)                                       :: uNum = 100 !Number of points in velocity space
    real(KREAL), parameter                              :: U_MIN = -15.0, U_MAX = +15.0 !Minimum and maximum micro velocity
    real(KREAL), allocatable, dimension(:)              :: uSpace !Discrete velocity space
    real(KREAL), allocatable, dimension(:)              :: weight !Qudrature weight at velocity u_k

end module ControlParameters

!--------------------------------------------------
!>Define some commonly used functions/subroutines
!--------------------------------------------------
module Tools
    use ControlParameters
    implicit none

contains
    !--------------------------------------------------
    !>Convert primary variables to conservative variables
    !>@param[in] prim          :primary variables
    !>@return    GetConserved  :conservative variables
    !--------------------------------------------------
    function GetConserved(prim)
        real(KREAL), intent(in)                         :: prim(3) !Density, x-velocity, lambda=1/temperature
        real(KREAL)                                     :: GetConserved(3) !Density, x-momentum, total energy

        GetConserved(1) = prim(1)
        GetConserved(2) = prim(1)*prim(2)
        GetConserved(3) = 0.5*prim(1)/(prim(3)*(GAMMA-1.0))+0.5*prim(1)*prim(2)**2
    end function GetConserved

    !--------------------------------------------------
    !>Convert conservative variables to primary variables
    !>@param[in] w           :conservative variables
    !>@return    GetPrimary  :primary variables
    !--------------------------------------------------
    function GetPrimary(w)
        real(KREAL), intent(in)                         :: w(3) !Density, x-momentum, total energy
        real(KREAL)                                     :: GetPrimary(3) !Density, x-velocity, lambda=1/temperature

        GetPrimary(1) = w(1)
        GetPrimary(2) = w(2)/w(1)
        GetPrimary(3) = 0.5*w(1)/((GAMMA-1.0)*(w(3)-0.5*w(2)**2/w(1)))
    end function GetPrimary

    !--------------------------------------------------
    !>Obtain speed of sound
    !>@param[in] prim    :primary variables
    !>@return    GetSoundSpeed :speed of sound
    !--------------------------------------------------
    function GetSoundSpeed(prim)
        real(KREAL), intent(in)                         :: prim(3)
        real(KREAL)                                     :: GetSoundSpeed !speed of sound

        GetSoundSpeed = sqrt(0.5*GAMMA/prim(3))
    end function GetSoundSpeed

    !--------------------------------------------------
    !>Obtain discretized Maxwellian distribution
    !>@param[out] h,b  :distribution function
    !>@param[in]  prim :primary variables
    !--------------------------------------------------
    subroutine DiscreteMaxwell(h,b,prim)
        real(KREAL), dimension(:), intent(out)          :: h,b
        real(KREAL), intent(in)                         :: prim(3)

        h = prim(1)*(prim(3)/PI)**(1.0/2.0)*exp(-prim(3)*(uSpace-prim(2))**2)
        b = h*CK/(2.0*prim(3))
    end subroutine DiscreteMaxwell

    !--------------------------------------------------
    !>Calculate the Shakhov part H^+, B^+
    !>@param[in]  H,B           :Maxwellian distribution function
    !>@param[in]  qf            :heat flux
    !>@param[in]  prim          :primary variables
    !>@param[out] H_plus,B_plus :Shakhov part
    !--------------------------------------------------
    subroutine ShakhovPart(H,B,qf,prim,H_plus,B_plus)
        real(KREAL), dimension(:), intent(in)           :: H,B
        real(KREAL), intent(in)                         :: qf
        real(KREAL), intent(in)                         :: prim(3)
        real(KREAL), dimension(:), intent(out)          :: H_plus,B_plus

        H_plus = 0.8*(1-PR)*prim(3)**2/prim(1)*&
                    (uSpace-prim(2))*qf*(2*prim(3)*(uSpace-prim(2))**2+CK-5)*H
        B_plus = 0.8*(1-PR)*prim(3)**2/prim(1)*&
                    (uSpace-prim(2))*qf*(2*prim(3)*(uSpace-prim(2))**2+CK-3)*B
    end subroutine ShakhovPart
    
    !--------------------------------------------------
    !>VanLeerLimiter for reconstruction of distrubution function
    !>@param[in]    leftCell  :the left cell
    !>@param[inout] midCell   :the middle cell
    !>@param[in]    rightCell :the right cell
    !--------------------------------------------------
    subroutine VanLeerLimiter(leftCell,midCell,rightCell)
        type(CellCenter), intent(in)                    :: leftCell,rightCell
        type(CellCenter), intent(inout)                 :: midCell
        real(KREAL), allocatable, dimension(:)          :: sL,sR

        !allocate array
        allocate(sL(uNum))
        allocate(sR(uNum))

        sL = (midCell%h-leftCell%h)/(0.5*(midCell%length+leftCell%length))
        sR = (rightCell%h-midCell%h)/(0.5*(rightCell%length+midCell%length))
        midCell%sh = (sign(UP,sR)+sign(UP,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)

        sL = (midCell%b-leftCell%b)/(0.5*(midCell%length+leftCell%length))
        sR = (rightCell%b-midCell%b)/(0.5*(rightCell%length+midCell%length))
        midCell%sb = (sign(UP,sR)+sign(UP,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)
    end subroutine VanLeerLimiter
end module Tools

!--------------------------------------------------
!>Flux calculation
!--------------------------------------------------
module Flux
    use Tools
    implicit none
    integer(KREAL), parameter                           :: MNUM = 6 !Number of normal velocity moments

contains
    !--------------------------------------------------
    !>Calculate flux of inner interface
    !>@param[in]    leftCell  :cell left to the target interface
    !>@param[inout] face      :the target interface
    !>@param[in]    rightCell :cell right to the target interface
    !--------------------------------------------------
    subroutine CalcFlux(leftCell,face,rightCell)
        type(CellCenter), intent(in)                    :: leftCell,rightCell
        type(CellInterface), intent(inout)              :: face
        real(KREAL), allocatable, dimension(:)          :: h,b !Distribution function at the interface
        real(KREAL), allocatable, dimension(:)          :: H0,B0 !Maxwellian distribution function
        real(KREAL), allocatable, dimension(:)          :: H_plus,B_plus !Shakhov part of the equilibrium distribution
        real(KREAL), allocatable, dimension(:)          :: sh,sb !Slope of distribution function at the interface
        integer(KINT), allocatable, dimension(:)        :: delta !Heaviside step function
        real(KREAL)                                     :: conVars(3),prim(3) !Conservative and primary variables at the interface
        real(KREAL)                                     :: qf !Heat flux in normal and tangential direction
        real(KREAL)                                     :: sw(3) !Slope of conVars
        real(KREAL)                                     :: aL(3),aR(3),aT(3) !Micro slope of Maxwellian distribution, left,right and time.
        real(KREAL)                                     :: Mu(0:MNUM),MuL(0:MNUM),MuR(0:MNUM),Mxi(0:2) !<u^n>,<u^n>_{>0},<u^n>_{<0},<\xi^l>
        real(KREAL)                                     :: Mau0(3),MauL(3),MauR(3),MauT(3) !<u\psi>,<aL*u^n*\psi>,<aR*u^n*\psi>,<A*u*\psi>
        real(KREAL)                                     :: tau !Collision time
        real(KREAL)                                     :: Mt(5) !Some time integration terms
        integer(KINT)                                   :: i,j

        !--------------------------------------------------
        !Prepare
        !--------------------------------------------------
        !Allocate array
        allocate(delta(uNum))
        allocate(h(uNum))
        allocate(b(uNum))
        allocate(sh(uNum))
        allocate(sb(uNum))
        allocate(H0(uNum))
        allocate(B0(uNum))
        allocate(H_plus(uNum))
        allocate(B_plus(uNum))

        !Heaviside step function
        delta = (sign(UP,uSpace)+1)/2

        !--------------------------------------------------
        !Reconstruct initial distribution at interface
        !--------------------------------------------------
        h = (leftCell%h+0.5*leftCell%length*leftCell%sh)*delta+&
            (rightCell%h-0.5*rightCell%length*rightCell%sh)*(1-delta)
        b = (leftCell%b+0.5*leftCell%length*leftCell%sb)*delta+&
            (rightCell%b-0.5*rightCell%length*rightCell%sb)*(1-delta)
        sh = leftCell%sh*delta+rightCell%sh*(1-delta)
        sb = leftCell%sb*delta+rightCell%sb*(1-delta)

        !--------------------------------------------------
        !Obtain macroscopic variables
        !--------------------------------------------------
        !Conservative variables w_0 at interface
        conVars(1) = sum(weight*h)
        conVars(2) = sum(weight*uSpace*h)
        conVars(3) = 0.5*(sum(weight*uSpace**2*h)+sum(weight*b))

        !Convert to primary variables
        prim = GetPrimary(conVars)

        !--------------------------------------------------
        !Calculate a^L,a^R
        !--------------------------------------------------
        sw = (conVars-leftCell%conVars)/(0.5*leftCell%length) !left slope of conVars
        aL = MicroSlope(prim,sw) !calculate a^L

        sw = (rightCell%conVars-conVars)/(0.5*rightCell%length) !right slope of conVars
        aR = MicroSlope(prim,sw) !calculate a^R

        !--------------------------------------------------
        !Calculate time slope of conVars and A
        !--------------------------------------------------
        !<u^n>,<\xi^l>,<u^n>_{>0},<u^n>_{<0}
        call CalcMoment(prim,Mu,Mxi,MuL,MuR) 

        MauL = Moment_au(aL,MuL,Mxi,1) !<aL*u*\psi>_{>0}
        MauR = Moment_au(aR,MuR,Mxi,1) !<aR*u*\psi>_{<0}

        sw = -prim(1)*(MauL+MauR) !Time slope of conVars
        aT = MicroSlope(prim,sw) !Calculate A

        !--------------------------------------------------
        !Calculate collision time and some time integration terms
        !--------------------------------------------------
        tau = GetTau(prim)

        Mt(4) = tau*(1.0-exp(-dt/tau))
        Mt(5) = -tau*dt*exp(-dt/tau)+tau*Mt(4)
        Mt(1) = dt-Mt(4)
        Mt(2) = -tau*Mt(1)+Mt(5) 
        Mt(3) = 0.5*dt**2-tau*Mt(1)

        !--------------------------------------------------
        !Calculate the flux of conservative variables related to g0
        !--------------------------------------------------
        Mau0 = Moment_uxi(Mu,Mxi,1,0) !<u*\psi>
        MauL = Moment_au(aL,MuL,Mxi,2) !<aL*u^2*\psi>_{>0}
        MauR = Moment_au(aR,MuR,Mxi,2) !<aR*u^2*\psi>_{<0}
        MauT = Moment_au(aT,Mu,Mxi,1) !<A*u*\psi>

        face%flux = Mt(1)*prim(1)*Mau0+Mt(2)*prim(1)*(MauL+MauR)+Mt(3)*prim(1)*MauT

        !--------------------------------------------------
        !calculate the flux of conservative variables related to g+ and f0
        !--------------------------------------------------
        !Maxwellian distribution H0 and B0
        call DiscreteMaxwell(H0,B0,prim)
    
        !Calculate heat flux
        qf = GetHeatFlux(h,b,prim) 

        !Shakhov part H+ and B+
        call ShakhovPart(H0,B0,qf,prim,H_plus,B_plus)

        !Conservative flux related to g+ and f0
        face%flux(1) = face%flux(1)+Mt(1)*sum(weight*uSpace*H_plus)+Mt(4)*sum(weight*uSpace*h)-Mt(5)*sum(weight*uSpace**2*sh)
        face%flux(2) = face%flux(2)+Mt(1)*sum(weight*uSpace**2*H_plus)+Mt(4)*sum(weight*uSpace**2*h)-Mt(5)*sum(weight*uSpace**3*sh)
        face%flux(3) = face%flux(3)+&
                        Mt(1)*0.5*(sum(weight*uSpace**3*H_plus)+sum(weight*uSpace*B_plus))+&
                        Mt(4)*0.5*(sum(weight*uSpace**3*h)+sum(weight*uSpace*b))-&
                        Mt(5)*0.5*(sum(weight*uSpace**4*sh)+sum(weight*uSpace**2*sb))

        !--------------------------------------------------
        !Calculate flux of distribution function
        !--------------------------------------------------
        face%flux_h = Mt(1)*uSpace*(H0+H_plus)+&
                        Mt(2)*uSpace**2*(aL(1)*H0+aL(2)*uSpace*H0+0.5*aL(3)*(uSpace**2*H0+B0))*delta+&
                        Mt(2)*uSpace**2*(aR(1)*H0+aR(2)*uSpace*H0+0.5*aR(3)*(uSpace**2*H0+B0))*(1-delta)+&
                        Mt(3)*uSpace*(aT(1)*H0+aT(2)*uSpace*H0+0.5*aT(3)*(uSpace**2*H0+B0))+&
                        Mt(4)*uSpace*h-Mt(5)*uSpace**2*sh

        face%flux_b = Mt(1)*uSpace*(B0+B_plus)+&
                        Mt(2)*uSpace**2*(aL(1)*B0+aL(2)*uSpace*B0+0.5*aL(3)*(uSpace**2*B0+Mxi(2)*H0))*delta+&
                        Mt(2)*uSpace**2*(aR(1)*B0+aR(2)*uSpace*B0+0.5*aR(3)*(uSpace**2*B0+Mxi(2)*H0))*(1-delta)+&
                        Mt(3)*uSpace*(aT(1)*B0+aT(2)*uSpace*B0+0.5*aT(3)*(uSpace**2*B0+Mxi(2)*H0))+&
                        Mt(4)*uSpace*b-Mt(5)*uSpace**2*sb
        
        !--------------------------------------------------
        !Aftermath
        !--------------------------------------------------
        !Deallocate array
        deallocate(delta)
        deallocate(h)
        deallocate(b)
        deallocate(sh)
        deallocate(sb)
        deallocate(H0)
        deallocate(B0)
        deallocate(H_plus)
        deallocate(B_plus)
    end subroutine CalcFlux

    !--------------------------------------------------
    !>Calculate micro slope of Maxwellian distribution
    !>@param[in] prim        :primary variables
    !>@param[in] sw          :slope of W
    !>@return    MicroSlope  :slope of Maxwellian distribution
    !--------------------------------------------------
    function MicroSlope(prim,sw)
        real(KREAL), intent(in)                         :: prim(3),sw(3)
        real(KREAL)                                     :: MicroSlope(3)

        MicroSlope(3) = 4.0*prim(3)**2/((CK+1)*prim(1))*(2.0*sw(3)-2.0*prim(2)*sw(2)+sw(1)*(prim(2)**2-0.5*(CK+1)/prim(3)))
        MicroSlope(2) = 2.0*prim(3)/prim(1)*(sw(2)-prim(2)*sw(1))-prim(2)*MicroSlope(3)
        MicroSlope(1) = sw(1)/prim(1)-prim(2)*MicroSlope(2)-0.5*(prim(2)**2+0.5*(CK+1)/prim(3))*MicroSlope(3)
    end function MicroSlope

    !--------------------------------------------------
    !>Calculate collision time
    !>@param[in] prim    :primary variables
    !>@return    GetTau  :collision time
    !--------------------------------------------------
    function GetTau(prim)
        real(KREAL), intent(in)                         :: prim(3)
        real(KREAL)                                     :: GetTau

        GetTau = MU_REF*2.0*prim(3)**(1-OMEGA)/prim(1)
    end function GetTau
    
    !--------------------------------------------------
    !>Calculate heat flux
    !>@param[in] h,b           :distribution function
    !>@param[in] prim          :primary variables
    !>@return    GetHeatFlux   :heat flux in normal and tangential direction
    !--------------------------------------------------
    function GetHeatFlux(h,b,prim)
        real(KREAL), dimension(:), intent(in)           :: h,b
        real(KREAL), intent(in)                         :: prim(3)
        real(KREAL)                                     :: GetHeatFlux !heat flux in normal and tangential direction

        GetHeatFlux = 0.5*(sum(weight*(uSpace-prim(2))*(uSpace-prim(2))**2*h)+sum(weight*(uSpace-prim(2))*b)) 
    end function GetHeatFlux

    !--------------------------------------------------
    !>calculate moments of velocity and \xi
    !>@param[in] prim :primary variables
    !>@param[out] Mu        :<u^n>
    !>@param[out] Mxi       :<\xi^2n>
    !>@param[out] MuL,MuR   :<u^n>_{>0},<u^n>_{<0}
    !--------------------------------------------------
    subroutine CalcMoment(prim,Mu,Mxi,MuL,MuR)
        real(KREAL), intent(in)                         :: prim(3)
        real(KREAL), intent(out)                        :: Mu(0:MNUM),MuL(0:MNUM),MuR(0:MNUM)
        real(KREAL), intent(out)                        :: Mxi(0:2)
        integer :: i

        !Moments of normal velocity
        MuL(0) = 0.5*erfc(-sqrt(prim(3))*prim(2))
        MuL(1) = prim(2)*MuL(0)+0.5*exp(-prim(3)*prim(2)**2)/sqrt(PI*prim(3))
        MuR(0) = 0.5*erfc(sqrt(prim(3))*prim(2))
        MuR(1) = prim(2)*MuR(0)-0.5*exp(-prim(3)*prim(2)**2)/sqrt(PI*prim(3))

        do i=2,MNUM
            MuL(i) = prim(2)*MuL(i-1)+0.5*(i-1)*MuL(i-2)/prim(3)
            MuR(i) = prim(2)*MuR(i-1)+0.5*(i-1)*MuR(i-2)/prim(3)
        end do

        Mu = MuL+MuR

        !Moments of \xi
        Mxi(0) = 1.0 !<\xi^0>
        Mxi(1) = 0.5*CK/prim(3) !<\xi^2>
        Mxi(2) = CK*(CK+2.0)/(4.0*prim(3)**2) !<\xi^4>
    end subroutine CalcMoment

    !--------------------------------------------------
    !>Calculate <a*u^\alpha*\psi>
    !>@param[in] a         :micro slope of Maxwellian
    !>@param[in] Mu        :<u^\alpha>
    !>@param[in] Mxi       :<\xi^l>
    !>@param[in] alpha     :exponential index of u
    !>@return    Moment_au :moment of <a*u^\alpha*\psi>
    !--------------------------------------------------
    function Moment_au(a,Mu,Mxi,alpha)
        real(KREAL), intent(in)                         :: a(3)
        real(KREAL), intent(in)                         :: Mu(0:MNUM),Mxi(0:2)
        integer(KINT), intent(in)                       :: alpha
        real(KREAL)                                     :: Moment_au(3)

        Moment_au = a(1)*Moment_uxi(Mu,Mxi,alpha+0,0)+&
                    a(2)*Moment_uxi(Mu,Mxi,alpha+1,0)+&
                    0.5*a(3)*Moment_uxi(Mu,Mxi,alpha+2,0)+&
                    0.5*a(3)*Moment_uxi(Mu,Mxi,alpha+0,2)
    end function Moment_au

    !--------------------------------------------------
    !>Calculate <u^\alpha*\xi^\delta*\psi>
    !>@param[in] Mu         :<u^\alpha>
    !>@param[in] Mxi        :<\xi^delta>
    !>@param[in] alpha      :exponential index of u
    !>@param[in] delta      :exponential index of \xi
    !>@return    Moment_uxi :moment of <u^\alpha*\xi^\delta*\psi>
    !--------------------------------------------------
    function Moment_uxi(Mu,Mxi,alpha,delta)
        real(KREAL), intent(in)                         :: Mu(0:MNUM),Mxi(0:2)
        integer(KINT), intent(in)                       :: alpha,delta
        real(KREAL)                                     :: Moment_uxi(3)

        Moment_uxi(1) = Mu(alpha)*Mxi(delta/2)
        Moment_uxi(2) = Mu(alpha+1)*Mxi(delta/2)
        Moment_uxi(3) = 0.5*(Mu(alpha+2)*Mxi(delta/2)+Mu(alpha)*Mxi((delta+2)/2))
    end function Moment_uxi
end module Flux

!--------------------------------------------------
!>UGKS solver
!--------------------------------------------------
module Solver
    use Flux
    implicit none

contains
    !--------------------------------------------------
    !>Calculate time step
    !--------------------------------------------------
    subroutine TimeStep()
        real(KREAL)                                     :: prim(3) !primary variables
        real(KREAL)                                     :: sos !speed of sound
        real(KREAL)                                     :: tMax
        integer(KINT)                                   :: i
        
        !Set initial value
        tMax = 0.0

        do i=IXMIN,IXMAX
            !Convert conservative variables to primary variables
            prim = GetPrimary(ctr(i)%conVars)

            !Get sound speed
            sos = GetSoundSpeed(prim)

            !Maximum velocity
            prim(2) = max(U_MAX,abs(prim(2)))+sos

            !Maximum 1/dt allowed
            tMax = max(tMax,prim(2)/ctr(i)%length)
        end do

        !Time step
        dt = CFL/tMax
    end subroutine TimeStep
    
    !--------------------------------------------------
    !Enforce boundary condition
    !--------------------------------------------------
    subroutine Boundary()
        !Get conservative variables and distribution function
        ctr(IXMIN-GHOST_NUM)%conVars = GetConserved(PRIM_LEFT)
        ctr(IXMAX+GHOST_NUM)%conVars = GetConserved(PRIM_RIGHT)
        call DiscreteMaxwell(ctr(IXMIN-GHOST_NUM)%h,ctr(IXMIN-GHOST_NUM)%b,PRIM_LEFT)
        call DiscreteMaxwell(ctr(IXMAX+GHOST_NUM)%h,ctr(IXMAX+GHOST_NUM)%b,PRIM_RIGHT)

        !Initialize slope of distribution function at ghost cell
        ctr(IXMIN-GHOST_NUM)%sh = 0.0; ctr(IXMIN-GHOST_NUM)%sb = 0.0
        ctr(IXMAX+GHOST_NUM)%sh = 0.0; ctr(IXMAX+GHOST_NUM)%sb = 0.0
    end subroutine Boundary

    !--------------------------------------------------
    !>Reconstruct the slope of initial distribution function
    !>Index method
    ! -------------------------------
    ! | (i-1) |  (i) |  (i) | (i+1) |
    ! | cell  | face | cell | face  |
    ! -------------------------------
    !--------------------------------------------------
    subroutine Reconstruction()
        integer(KINT) :: i

        do i=IXMIN-GHOST_NUM+1,IXMAX+GHOST_NUM-1
            call VanLeerLimiter(ctr(i-1),ctr(i),ctr(i+1))
        end do
    end subroutine Reconstruction

    !--------------------------------------------------
    !>Calculate the flux across the interfaces
    !>Index method
    ! -------------------------------
    ! | (i-1) |  (i) |  (i) | (i+1) |
    ! | cell  | face | cell | face  |
    ! -------------------------------
    !--------------------------------------------------
    subroutine Evolution()
        integer(KINT) :: i

        do i=IXMIN,IXMAX+1
            call CalcFlux(ctr(i-1),vface(i),ctr(i))
        end do
    end subroutine Evolution

    !--------------------------------------------------
    !>Update cell averaged values
    !>Index method
    ! -------------------------------
    ! | (i-1) |  (i) |  (i) | (i+1) |
    ! | cell  | face | cell | face  |
    ! -------------------------------
    !--------------------------------------------------
    subroutine Update()
        real(KREAL), allocatable, dimension(:)          :: H_old,B_old !Equilibrium distribution at t=t^n
        real(KREAL), allocatable, dimension(:)          :: H,B !Equilibrium distribution at t=t^{n+1}
        real(KREAL), allocatable, dimension(:)          :: H_plus,B_plus !Shakhov part
        real(KREAL)                                     :: prim_old(3),prim(3) !Primary variables at t^n and t^{n+1}
        real(KREAL)                                     :: tau_old,tau !Collision time and t^n and t^{n+1}
        real(KREAL)                                     :: qf
        integer(KINT)                                   :: i

        !Allocate arrays
        allocate(H_old(uNum))
        allocate(B_old(uNum))
        allocate(H(uNum))
        allocate(B(uNum))
        allocate(H_plus(uNum))
        allocate(B_plus(uNum))

        do i=IXMIN,IXMAX
            !--------------------------------------------------
            !Store conVars^n and calculate H^n,B^n,\tau^n
            !--------------------------------------------------
            prim_old = GetPrimary(ctr(i)%conVars) !Convert to primary variables
            call DiscreteMaxwell(H_old,B_old,prim_old) !Calculate Maxwellian
            tau_old = GetTau(prim_old) !Calculate collision time \tau^n

            !--------------------------------------------------
            !Update conVars^{n+1} and Calculate H^{n+1},B^{n+1},\tau^{n+1}
            !--------------------------------------------------
            ctr(i)%conVars = ctr(i)%conVars+(vface(i)%flux-vface(i+1)%flux)/ctr(i)%length !Update conVars^{n+1}

            prim = GetPrimary(ctr(i)%conVars)
            call DiscreteMaxwell(H,B,prim)
            tau = GetTau(prim)

            !--------------------------------------------------
            !Shakhov part
            !--------------------------------------------------
            !Calculate heat flux at t=t^n
            qf = GetHeatFlux(ctr(i)%h,ctr(i)%b,prim_old) 

            !h^+ = H+H^+ at t=t^n
            call ShakhovPart(H_old,B_old,qf,prim_old,H_plus,B_plus) !H^+ and B^+
            H_old = H_old+H_plus !h^+
            B_old = B_old+B_plus !b^+

            !h^+ = H+H^+ at t=t^{n+1}
            call ShakhovPart(H,B,qf,prim,H_plus,B_plus)
            H = H+H_plus
            B = B+B_plus

            !--------------------------------------------------
            !Update distribution function
            !--------------------------------------------------
            ctr(i)%h = (ctr(i)%h+(vface(i)%flux_h-vface(i+1)%flux_h)/ctr(i)%length+&
                                0.5*dt*(H/tau+(H_old-ctr(i)%h)/tau_old))/(1.0+0.5*dt/tau)
            ctr(i)%b = (ctr(i)%b+(vface(i)%flux_b-vface(i+1)%flux_b)/ctr(i)%length+&
                                0.5*dt*(B/tau+(B_old-ctr(i)%b)/tau_old))/(1.0+0.5*dt/tau)
        end do

        !Deallocate arrays
        deallocate(H_old)
        deallocate(B_old)
        deallocate(H)
        deallocate(B)
        deallocate(H_plus)
        deallocate(B_plus)
    end subroutine Update
end module Solver

!--------------------------------------------------
!>Initialization of mesh and intial flow field
!--------------------------------------------------
module Initialization
    use ControlParameters
    use Tools
    implicit none

contains
    !--------------------------------------------------
    !>Main initialization subroutine
    !--------------------------------------------------
    subroutine Init()  
        ! call InitUniformMesh() !Initialize Uniform mesh
        call InitRandomMesh()
        call InitVelocityNewton(uNum) !Initialize discrete velocity space
        call InitFlowField(PRIM_LEFT,PRIM_RIGHT) !Set the initial value
    end subroutine Init
    
    !--------------------------------------------------
    !>Initialize uniform mesh
    !--------------------------------------------------
    subroutine InitUniformMesh()
        real(KREAL)                                     :: dx
        integer(KINT)                                   :: i

        !Cell length
        dx = (END_POINT-START_POINT)/(IXMAX-IXMIN+1)

        !Cell center (with ghost cell)
        forall(i=IXMIN-GHOST_NUM:IXMAX+GHOST_NUM) 
            ctr(i)%x = START_POINT+(i-IXMIN+0.5)*dx
            ctr(i)%length = dx
        end forall
    end subroutine InitUniformMesh

    !--------------------------------------------------
    !>Initialize random mesh
    !--------------------------------------------------
    subroutine InitRandomMesh()
        real(KREAL)                                     :: dx(IXMIN-GHOST_NUM:IXMAX+GHOST_NUM)
        integer(KINT)                                   :: i

        call random_number(dx)

        !Cell length
        dx(IXMIN:IXMAX) = dx(IXMIN:IXMAX)/sum(dx(IXMIN:IXMAX))*(END_POINT-START_POINT)
        dx(IXMIN-GHOST_NUM:IXMIN) = dx(IXMIN); dx(IXMAX:IXMAX+GHOST_NUM) = dx(IXMAX)

        !Cell center (with ghost cell)
        ctr(IXMIN-GHOST_NUM)%x = START_POINT-dx(IXMIN)*GHOST_NUM; ctr(IXMIN-GHOST_NUM)%length = dx(IXMIN-GHOST_NUM)
        do i=IXMIN-GHOST_NUM+1,IXMAX+GHOST_NUM
            ctr(i)%x = ctr(i-1)%x+dx(i)
            ctr(i)%length = dx(i)
        end do
    end subroutine InitRandomMesh

    !--------------------------------------------------
    !>Initialize discrete velocity space using Newton–Cotes formulas
    !--------------------------------------------------
    subroutine InitVelocityNewton(num)
        integer(KINT), intent(inout)                    :: num
        real(KREAL)                                     :: du !Spacing in u velocity space
        integer(KINT)                                   :: i

        !Modify num if not appropriate
        num = (num/4)*4+1

        !Allocate array
        allocate(uSpace(num))
        allocate(weight(num))

        !spacing in u and v velocity space
        du = (U_MAX-U_MIN)/(num-1)

        !velocity space
        forall(i=1:num)
            uSpace(i) = U_MIN+(i-1)*du
            weight(i) = (newtonCoeff(i,num)*du)
        end forall

        !--------------------------------------------------
        !>Allocate arrays in velocity space
        !--------------------------------------------------
        !cell center (with ghost cell)
        do i=IXMIN-GHOST_NUM,IXMAX+GHOST_NUM
            allocate(ctr(i)%h(num))
            allocate(ctr(i)%b(num))
            allocate(ctr(i)%sh(num))
            allocate(ctr(i)%sb(num))
        end do

        !cell interface
        do i=IXMIN-GHOST_NUM+1,IXMAX+GHOST_NUM
            allocate(vface(i)%flux_h(num))
            allocate(vface(i)%flux_b(num))
        end do

        contains
            !--------------------------------------------------
            !>Calculate the coefficient for newton-cotes formula
            !>@param[in] idx          :index in velocity space
            !>@param[in] num          :total number in velocity space
            !>@return    newtonCoeff  :coefficient for newton-cotes formula
            !--------------------------------------------------
            pure function newtonCoeff(idx,num)
                integer(KINT), intent(in)               :: idx,num
                real(KREAL)                             :: newtonCoeff

                if (idx==1 .or. idx==num) then 
                    newtonCoeff = 14.0/45.0
                else if (mod(idx-5,4)==0) then
                    newtonCoeff = 28.0/45.0
                else if (mod(idx-3,4)==0) then
                    newtonCoeff = 24.0/45.0
                else
                    newtonCoeff = 64.0/45.0
                end if
            end function newtonCoeff
    end subroutine InitVelocityNewton

    !--------------------------------------------------
    !>Set the initial condition
    !--------------------------------------------------
    subroutine InitFlowField(primL,primR)
        real(KREAL), intent(in)                         :: primL(3), primR(3) !primary variables before and after shock
        real(KREAL)                                     :: wL(3), wR(3) !conservative variables before and after shock
        real(KREAL), allocatable, dimension(:)          :: hL,bL !distribution function before shock
        real(KREAL), allocatable, dimension(:)          :: hR,bR !distribution function after shock
        integer(KINT)                                   :: i

        !Allocation
        allocate(hL(uNum))
        allocate(bL(uNum))
        allocate(hR(uNum))
        allocate(bR(uNum))

        !Get conservative variables and distribution function
        wL = GetConserved(primL)
        wR = GetConserved(primR)
        call DiscreteMaxwell(hL,bL,primL)
        call DiscreteMaxwell(hR,bR,primR)

        !Initialize field (with ghost cell)
        forall(i=IXMIN-GHOST_NUM:(IXMIN+IXMAX)/2)
            ctr(i)%conVars = wL
            ctr(i)%h = hL
            ctr(i)%b = bL
        end forall
        
        forall(i=(IXMIN+IXMAX)/2+1:IXMAX+GHOST_NUM)
            ctr(i)%conVars = wR
            ctr(i)%h = hR
            ctr(i)%b = bR
        end forall 

        !Initialize slope of distribution function at ghost cell
        ctr(IXMIN-GHOST_NUM)%sh = 0.0; ctr(IXMIN-GHOST_NUM)%sb = 0.0
        ctr(IXMAX+GHOST_NUM)%sh = 0.0; ctr(IXMAX+GHOST_NUM)%sb = 0.0

        !Deallocation
        deallocate(hL)
        deallocate(bL)
        deallocate(hR)
        deallocate(bR)
    end subroutine InitFlowField
    
end module Initialization

module Writer
    use Tools
    implicit none
    character(len=8)                                    :: date
    character(len=10)                                   :: time
    character(len=100)                                  :: fileName
contains
    !--------------------------------------------------
    !>Write result
    !--------------------------------------------------
    subroutine Output()
        integer(KINT)                                   :: i
        real(KREAL)                                     :: rmid !Average density
        real(KREAL)                                     :: xmid !Location of average density
        real(KREAL), dimension(:,:), allocatable        :: solution
        !--------------------------------------------------
        !Prepare solutions
        !--------------------------------------------------
        allocate(solution(3,IXMIN-GHOST_NUM:IXMAX+GHOST_NUM)) !including ghost cell

        do i=IXMIN-GHOST_NUM,IXMAX+GHOST_NUM
            solution(:,i) = GetPrimary(ctr(i)%conVars)
            solution(3,i) = 1.0/solution(3,i) !Temperature=1/Lambda
        end do
        
        !Find middle location - the location of average density
        rmid = 0.5*(solution(1,IXMIN-1)+solution(1,IXMAX+1))

        do i=IXMIN,IXMAX
            if ((solution(1,i)-rmid)*(solution(1,i+1)-rmid)<=0) then
                xmid = ctr(i)%x+(ctr(i+1)%x-ctr(i)%x)/(solution(1,i+1)-solution(1,i))*(rmid-solution(1,i))
            end if
        end do

        !Normalization
        do i=1,3
            solution(i,:) = (solution(i,:)-solution(i,IXMIN-GHOST_NUM))/(solution(i,IXMAX+GHOST_NUM)-solution(i,IXMIN-GHOST_NUM))
        end do

        !--------------------------------------------------
        !Write to file
        !--------------------------------------------------
        !Open result file and write header
        ! using keyword arguments
        open(unit=RSTFILE,file=RSTFILENAME//trim(fileName)//'.dat',status="replace",action="write")
        ! write(RSTFILE,*) "VARIABLES = x, Density, U, Temperature"
        write(RSTFILE,*) "VARIABLES = x, Density, U, Temperature"
        write(RSTFILE,*) 'ZONE  T="Time: ',simTime,'", I = ',IXMAX-IXMIN+1,', DATAPACKING=BLOCK'

        !Write geometry (cell-centered)
        write(RSTFILE,"(6(ES23.16,2X))") ctr(IXMIN:IXMAX)%x-xmid

        !Write solution (cell-centered)
        do i=1,3
            write(RSTFILE,"(6(ES23.16,2X))") solution(i,IXMIN:IXMAX)
        end do

        !close file
        close(RSTFILE)
        deallocate(solution)
    end subroutine Output
end module Writer
!--------------------------------------------------
!>main program
!--------------------------------------------------
program StationaryShockStructure
    use Initialization
    use Solver
    use Writer
    implicit none
    real(KREAL)                                         :: start, finish
    integer(KINT)                                       :: i
    
    !Initialization
    call Init()

    !Open file and write header
    call date_and_time(DATE=date,TIME=time)
    fileName = '_'//date//'_'//time(1:6)
    open(unit=HSTFILE,file=HSTFILENAME//trim(fileName)//'.hst',status="replace",action="write") !open history file
    write(HSTFILE,*) "VARIABLES = iter, simTime, dt" !write header

    !Star timer
    call cpu_time(start)

    !Iteration
    do while(.true.)
        call TimeStep() !Calculate the time step
        call Boundary() !Set up boundary condition
        call Reconstruction() !Initial value reconstruction
        call Evolution() !Calculate flux across the interfaces
        call Update() !Update cell averaged value

        !Check stopping criterion
        if (simTime>=MAX_TIME .or. iter>=MAX_ITER) exit

        !Log the iteration situation every 10 iterations
        if (mod(iter,10)==0) then
            write(*,"(A18,I15,2E15.7)") "iter,simTime,dt:",iter,simTime,dt
            write(HSTFILE,"(I15,2E15.7)") iter,simTime,dt
        end if

        iter = iter+1
        simTime = simTime+dt
    enddo

    !End timer
    call cpu_time(finish)
    print '("Run Time = ",f6.3," seconds.")', finish-start

    !Close history file
    close(HSTFILE)

    !Output solution
    call Output()

    !Aftermath
    !Deallocate array
    deallocate(uSpace)
    deallocate(weight)
    !Cell center
    do i=IXMIN-GHOST_NUM,IXMAX+GHOST_NUM
        deallocate(ctr(i)%h)
        deallocate(ctr(i)%b)
        deallocate(ctr(i)%sh)
        deallocate(ctr(i)%sb)
    end do
    !Cell interface
    do i=IXMIN-GHOST_NUM+1,IXMAX+GHOST_NUM
        deallocate(vface(i)%flux_h)
        deallocate(vface(i)%flux_b)
    end do
end program StationaryShockStructure