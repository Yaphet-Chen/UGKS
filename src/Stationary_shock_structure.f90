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
    subroutine VanleerLimiter(leftCell,midCell,rightCell,midFace)
        type(CellCenter), intent(in)                    :: leftCell, rightCell
        type(CellCenter), intent(inout)                 :: midCell
        type(CellInterface), intent(inout)              :: midface
        real(KREAL), dimension(3)                       :: splus, sminus
        real(KREAL), parameter                          :: UP = 1.0_KREAL

        splus = (rightCell%conVars-midCell%conVars)/(0.5*(rightCell%length+midCell%length))
        sminus = (midCell%conVars-leftCell%conVars)/(0.5*(midCell%length+leftCell%length))

        midface%leftVars = midCell%conVars-0.5*midCell%length*(sign(UP,splus)+sign(UP,sminus))*abs(splus)*abs(sminus)/(abs(splus)+abs(sminus)+TINYNUM)
        midface%rightVars = midCell%conVars+0.5*midCell%length*(sign(UP,splus)+sign(UP,sminus))*abs(splus)*abs(sminus)/(abs(splus)+abs(sminus)+TINYNUM)
        
        ! avoid blow up
        if ( GetLambda(midCell%leftVars)<=SMV .or. GetLambda(midCell%rightVars)<=SMV ) then
            midCell%leftVars = midCell%conVars
            midCell%rightVars = midCell%conVars
        endif
    end subroutine VanleerLimiter
end module Tools
!--------------------------------------------------
!>flux calculation
!--------------------------------------------------
module Flux
    use Tools
    implicit none

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
program Stationary_shock_structure
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

end program Stationary_shock_structure