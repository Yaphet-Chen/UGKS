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

    real(KREAL), parameter                              :: PI = 4.0_KREAL*atan(1.0_KREAL) !Pi
    real(KREAL), parameter                              :: SMV = tiny(0.0_KREAL) !Small value to avoid 0/0
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
        !Flow field
        real(KREAL)                                     :: leftVars(3) !Left conservative variables at cell interface
        real(KREAL)                                     :: rightVars(3) !Right conservative variables at cell interface
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
    real(KREAL), parameter                              :: CFL = 0.95_KREAL !CFL number
    real(KREAL)                                         :: simTime = 0.0_KREAL !Current simulation time
    real(KREAL), parameter                              :: MAX_TIME = 250.0_KREAL !Maximal simulation time
    integer(KINT), parameter                            :: MAX_ITER = 100 !Maximal iteration number
    integer(KINT)                                       :: iter = 1_KINT !Number of iteration
    real(KREAL)                                         :: dt !Global time step

    !Output control
    character(len=28), parameter                        :: HSTFILENAME = "StationaryShockStructure.hst" !History file name
    character(len=28), parameter                        :: RSTFILENAME = "StationaryShockStructure.dat" !Result file name
    integer, parameter                                  :: HSTFILE = 20 !History file ID
    integer, parameter                                  :: RSTFILE = 21 !Result file ID

    !Gas propeties
    integer(KINT), parameter                            :: CK = 2 !Internal degree of freedom, here 2 denotes monatomic gas
    real(KREAL), parameter                              :: GAMMA = real(k+3,KREAL)/real(k+1,KREAL) !Ratio of specific heat
    real(KREAL), parameter                              :: OMEGA = 0.72 !Temperature dependence index in HS/VHS/VSS model
    real(KREAL), parameter                              :: PR = 2.0/3.0 !Prandtl number
    real(KREAL), parameter                              :: KN = 1.0 !Knudsen number in reference state
    real(KREAL), parameter                              :: ALPHA_REF = 1.0 !Coefficient in HS model
    real(KREAL), parameter                              :: OMEGA_REF = 0.5 !Coefficient in HS model
    real(KREAL), parameter                              :: MU_REF !Viscosity coefficient in reference state
    real(KREAL), parameter                              :: MA = 8.0 !Mach number
    !Geometry
    real(KREAL), parameter                              :: START_POINT = 0.0_KREAL, END_POINT = 50.0_KREAL
    integer(KINT), parameter                            :: POINTS_NUM = 100_KINT
    integer(KINT), parameter                            :: IXMIN = 1 , IXMAX = POINTS_NUM !Cell index range
    integer(KINT), parameter                            :: GHOST_NUM = 1 !Ghost cell number

    !--------------------------------------------------
    !Initial flow field
    !--------------------------------------------------
    !Index method
    ! --------------------
    !  (i) |  (i)  | (i+1)
    ! face | cell  | face
    ! --------------------
    type(CellCenter)                                    :: ctr(IXMIN-GHOST_NUM:IXMAX+GHOST_NUM) !Cell center (with ghost cell)
    type(CellInterface)                                 :: vface(IXMIN-GHOST_NUM+1:IXMAX+GHOST_NUM) !Vertical cell interfaces

    !Initial value
    !Upstream condition (before shock)
    real(KREAL), parameter                              :: PRIM_LEFT(1) = 1.0 !Density
    real(KREAL), parameter                              :: PRIM_LEFT(2) = MA*sqrt(GAMMA/2.0) !X-velocity
    real(KREAL), parameter                              :: PRIM_LEFT(3) = 1.0 !Lambda=1/temperature
    !Downstream condition (after shock)
    real(KREAL), parameter                              :: MA_R = sqrt((MA**2*(GAMMA-1)+2)/(2*GAMMA*MA**2-(GAMMA-1)))
    real(KREAL), parameter                              :: RATIO_T = (1+(GAMMA-1)/2*MA**2)*(2*GAMMA/(GAMMA-1)*MA**2-1)/(MA**2*(2*GAMMA/(GAMMA-1)+(GAMMA-1)/2))
    real(KREAL), parameter                              :: PRIM_RIGHT(1) = PRIM_LEFT(1)*(GAMMA+1)*MA**2/((GAMMA-1)*MA**2+2)
    real(KREAL), parameter                              :: PRIM_RIGHT(2) = MA_R*sqrt(GAMMA/2.0)*sqrt(RATIO_T)
    real(KREAL), parameter                              :: PRIM_RIGHT(3) = PRIM_LEFT(3)/RATIO_T

    !--------------------------------------------------
    !Discrete velocity space
    !--------------------------------------------------
    integer(KINT)                                       :: uNum = 100 !Number of points in velocity space
    real(KREAL), parameter                              :: U_MIN = -15.0, U_MAX = 15.0 !Minimum and maximum micro velocity
    real(KREAL), allocatable, dimension(:)              :: uSpace !discrete velocity space
    real(KREAL), allocatable, dimension(:)              :: weight !weight at velocity u_k
end module ControlParameters

!--------------------------------------------------
!>Define some commonly used functions/subroutines
!--------------------------------------------------
module Tools
    use ConstantVariables
    use Mesh
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
    !>Obtain discretized Maxwellian distribution
    !>@param[out] h,b  :distribution function
    !>@param[in]  prim :primary variables
    !--------------------------------------------------
    subroutine DiscreteMaxwell(h,b,prim)
        real(KREAL), dimension(:), intent(out)          :: h,b
        real(KREAL), intent(in)                         :: prim(3)

        h = prim(1)*(prim(3)/PI)**(1.0/2.0)*exp(-prim(3)*(uspace-prim(2))**2)
        b = h*CK/(2.0*prim(3))
    end subroutine DiscreteMaxwell

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
!>flux calculation
!--------------------------------------------------
module Flux
    use Tools
    implicit none
contains
    !--------------------------------------------------
    !>Calculate flux of inner interface
    !>@param[in]    leftCell  :cell left to the target interface
    !>@param[inout] face      :the target interface
    !>@param[in]    rightCell :cell right to the target interface
    !--------------------------------------------------
    subroutine CalcFlux(leftCell,face,rightCell)
        type(CellCenter), intent(in) :: leftCell,rightCell
        type(CellInterface), intent(inout) :: face
        real(KREAL), allocatable, dimension(:) :: h,b !distribution function at the interface
        real(KREAL), allocatable, dimension(:) :: H0,B0 !Maxwellian distribution function
        real(KREAL), allocatable, dimension(:) :: H_plus,B_plus !Shakhov part of the equilibrium distribution
        real(KREAL), allocatable, dimension(:) :: sh,sb !slope of distribution function at the interface
        integer(KINT), allocatable, dimension(:) :: delta !Heaviside step function
        real(KREAL) :: w(3),prim(3) !conservative and primary variables at the interface
        real(KREAL) :: qf !heat flux in normal and tangential direction
        real(KREAL) :: sw(3) !slope of W
        real(KREAL) :: aL(3),aR(3),aT(3) !micro slope of Maxwellian distribution, left,right and time.
        real(KREAL) :: Mu(0:MNUM),Mu_L(0:MNUM),Mu_R(0:MNUM),Mxi(0:2) !<u^n>,<u^n>_{>0},<u^n>_{<0},<\xi^l>
        real(KREAL) :: Mau_0(3),Mau_L(3),Mau_R(3),Mau_T(3) !<u\psi>,<aL*u^n*\psi>,<aR*u^n*\psi>,<A*u*\psi>
        real(KREAL) :: tau !collision time
        real(KREAL) :: Mt(5) !some time integration terms
        integer :: i,j

        !--------------------------------------------------
        !prepare
        !--------------------------------------------------
        !allocate array
        allocate(delta(unum))
        allocate(h(unum))
        allocate(b(unum))
        allocate(sh(unum))
        allocate(sb(unum))
        allocate(H0(unum))
        allocate(B0(unum))
        allocate(H_plus(unum))
        allocate(B_plus(unum))

        !Heaviside step function
        delta = (sign(UP,uspace)+1)/2

        !--------------------------------------------------
        !reconstruct initial distribution
        !--------------------------------------------------
        h = (leftCell%h+0.5*leftCell%length*leftCell%sh)*delta+&
            (rightCell%h-0.5*rightCell%length*rightCell%sh)*(1-delta)
        b = (leftCell%b+0.5*leftCell%length*leftCell%sb)*delta+&
            (rightCell%b-0.5*rightCell%length*rightCell%sb)*(1-delta)
        sh = leftCell%sh*delta+rightCell%sh*(1-delta)
        sb = leftCell%sb*delta+rightCell%sb*(1-delta)

        !--------------------------------------------------
        !obtain macroscopic variables
        !--------------------------------------------------
        !conservative variables w_0
        w(1) = sum(weight*h)
        w(2) = sum(weight*uspace*h)
        w(3) = 0.5*(sum(weight*uspace**2*h)+sum(weight*b))

        !convert to primary variables
        prim = get_primary(w)

        !heat flux
        qf = get_heat_flux(h,b,prim) 

        !--------------------------------------------------
        !calculate a^L,a^R
        !--------------------------------------------------
        sw = (w-leftCell%w)/(0.5*leftCell%length) !left slope of W
        aL = micro_slope(prim,sw) !calculate a^L

        sw = (rightCell%w-w)/(0.5*rightCell%length) !right slope of W
        aR = micro_slope(prim,sw) !calculate a^R

        !--------------------------------------------------
        !calculate time slope of W and A
        !--------------------------------------------------
        !<u^n>,<\xi^l>,<u^n>_{>0},<u^n>_{<0}
        call calc_moment_u(prim,Mu,Mxi,Mu_L,Mu_R) 

        Mau_L = moment_au(aL,Mu_L,Mxi,1) !<aL*u*\psi>_{>0}
        Mau_R = moment_au(aR,Mu_R,Mxi,1) !<aR*u*\psi>_{<0}

        sw = -prim(1)*(Mau_L+Mau_R) !time slope of W
        aT = micro_slope(prim,sw) !calculate A

        !--------------------------------------------------
        !calculate collision time and some time integration terms
        !--------------------------------------------------
        tau = get_tau(prim)

        Mt(4) = tau*(1.0-exp(-dt/tau))
        Mt(5) = -tau*dt*exp(-dt/tau)+tau*Mt(4)
        Mt(1) = dt-Mt(4)
        Mt(2) = -tau*Mt(1)+Mt(5) 
        Mt(3) = dt**2/2.0-tau*Mt(1)

        !--------------------------------------------------
        !calculate the flux of conservative variables related to g0
        !--------------------------------------------------
        Mau_0 = moment_uv(Mu,Mxi,1,0) !<u*\psi>
        Mau_L = moment_au(aL,Mu_L,Mxi,2) !<aL*u^2*\psi>_{>0}
        Mau_R = moment_au(aR,Mu_R,Mxi,2) !<aR*u^2*\psi>_{<0}
        Mau_T = moment_au(aT,Mu,Mxi,1) !<A*u*\psi>

        face%flux = Mt(1)*prim(1)*Mau_0+Mt(2)*prim(1)*(Mau_L+Mau_R)+Mt(3)*prim(1)*Mau_T

        !--------------------------------------------------
        !calculate the flux of conservative variables related to g+ and f0
        !--------------------------------------------------
        !Maxwellian distribution H0 and B0
        call discrete_maxwell(H0,B0,prim)

        !Shakhov part H+ and B+
        call shakhov_part(H0,B0,qf,prim,H_plus,B_plus)

        !macro flux related to g+ and f0
        face%flux(1) = face%flux(1)+Mt(1)*sum(weight*uspace*H_plus)+Mt(4)*sum(weight*uspace*h)-Mt(5)*sum(weight*uspace**2*sh)
        face%flux(2) = face%flux(2)+Mt(1)*sum(weight*uspace**2*H_plus)+Mt(4)*sum(weight*uspace**2*h)-Mt(5)*sum(weight*uspace**3*sh)
        face%flux(3) = face%flux(3)+&
                        Mt(1)*0.5*(sum(weight*uspace*uspace**2*H_plus)+sum(weight*uspace*B_plus))+&
                        Mt(4)*0.5*(sum(weight*uspace*uspace**2*h)+sum(weight*uspace*b))-&
                        Mt(5)*0.5*(sum(weight*uspace**2*uspace**2*sh)+sum(weight*uspace**2*sb))

        !--------------------------------------------------
        !calculate flux of distribution function
        !--------------------------------------------------
        face%flux_h = Mt(1)*uspace*(H0+H_plus)+&
                        Mt(2)*uspace**2*(aL(1)*H0+aL(2)*uspace*H0+0.5*aL(3)*(uspace**2*H0+B0))*delta+&
                        Mt(2)*uspace**2*(aR(1)*H0+aR(2)*uspace*H0+0.5*aR(3)*(uspace**2*H0+B0))*(1-delta)+&
                        Mt(3)*uspace*(aT(1)*H0+aT(2)*uspace*H0+0.5*aT(3)*(uspace**2*H0+B0))+&
                        Mt(4)*uspace*h-Mt(5)*uspace**2*sh

        face%flux_b = Mt(1)*uspace*(B0+B_plus)+&
                        Mt(2)*uspace**2*(aL(1)*B0+aL(2)*uspace*B0+0.5*aL(3)*(uspace**2*B0+Mxi(2)*H0))*delta+&
                        Mt(2)*uspace**2*(aR(1)*B0+aR(2)*uspace*B0+0.5*aR(3)*(uspace**2*B0+Mxi(2)*H0))*(1-delta)+&
                        Mt(3)*uspace*(aT(1)*B0+aT(2)*uspace*B0+0.5*aT(3)*(uspace**2*B0+Mxi(2)*H0))+&
                        Mt(4)*uspace*b-Mt(5)*uspace**2*sb
    end subroutine CalcFlux
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
        integer                                         :: i
        
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
    !>Reconstruct the slope of initial distribution function
    !--------------------------------------------------
    subroutine Reconstruction()
        integer(KINT) :: i

        do i=IXMIN-GHOST_NUM+1,IXMAX+GHOST_NUM-1
            call VanLeerLimiter(ctr(i-1),ctr(i),ctr(i+1))
        end do
    end subroutine Reconstruction

    !--------------------------------------------------
    !>Calculate the flux across the interfaces
    !--------------------------------------------------
    subroutine Evolution()
        integer(KINT) :: i

        do i=IXMIN,IXMAX+1
            call CalcFlux(ctr(i-1),vface(i),ctr(i))
        end do
    end subroutine Evolution
end module Solver
!--------------------------------------------------
!>Initialization of mesh and intial flow field
!--------------------------------------------------
module Initialization
    use ControlParameters
    implicit none

contains
    !--------------------------------------------------
    !>Main initialization subroutine
    !--------------------------------------------------
    subroutine Init()  
        call InitUniformMesh() !Initialize Uniform mesh
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
    !>Initialize discrete velocity space using Newton–Cotes formulas
    !--------------------------------------------------
    subroutine InitVelocityNewton(num)
        integer, intent(inout) :: num
        real(kind=RKD) :: du !Spacing in u velocity space
        integer :: i

        !modify num if not appropriate
        num = (num/4)*4+1

        !allocate array
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
        ctr(IXMIN-GHOST_NUM)%sh = 0.0
        ctr(IXMIN-GHOST_NUM)%sb = 0.0
        ctr(IXMAX+GHOST_NUM)%sh = 0.0
        ctr(IXMAX+GHOST_NUM)%sb = 0.0
    end subroutine InitFlowField
    
end module Initialization

!--------------------------------------------------
!>main program
!--------------------------------------------------
program StationaryShockStructure
    use Solver
    use Writer
    implicit none
    real(KREAL) :: start, finish
    
    !Initialization
    call Init()

    !Open file and write header
    open(unit=HSTFILE,file=HSTFILENAME,status="replace",action="write") !open history file
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

    !close history file
    close(HSTFILE)

    !output solution
    call output()

end program StationaryShockStructure