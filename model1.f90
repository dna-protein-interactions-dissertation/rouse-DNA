program model1
	
use randgen
implicit none

integer, parameter :: long = selected_real_kind(15, 50)
double precision, parameter :: PI=3.141592653589793238462
real(kind=long), parameter :: Dp = 2e-17
real(kind=long), parameter :: db = 0.01e-6
real(kind=long), parameter :: di = 1e-6
real(kind=long), parameter, dimension(10) :: d0= [0.01e-6, 0.02154e-6, 0.0464e-6, 0.1e-6, 0.2154e-6, &
	0.464e-6, 1e-6, 2.154e-6, 4.64e-6, 10e-6] 
integer, parameter :: N = 101
integer, parameter :: reps = 1000 ! number of times to run the test
real(kind=long), parameter :: dt = 0.8e-6
real(kind=long), parameter :: b = 60e-9
real(kind=long), parameter :: sigma = 1.2e-9
real(kind=long), parameter :: kB = 1.4e-23
real(kind=long), parameter :: T = 300
real(kind=long), parameter :: vs = 0.001
real(kind=long), parameter :: z = 6*PI*sigma*vs
real(kind=long), parameter :: k = (3*kB*T)/(b**2)

real(kind=long), dimension(3) :: sphere, p
real(kind=long), dimension(2) :: seed
real(kind=long), dimension(N, 3) :: r, rnew
real(kind=long), dimension(N+2, 3) :: s
real(kind=long) :: min_distance, phi, lambda, suml=0
logical :: sts

integer :: i, j, x, tot
	
	
do y = 1,10:		
	tot = 0

	do x = 1, reps

		sts = .TRUE. ! This is the status of the simulation - should we keep going or not
		
		! We begin by initializing the DNA via the FJC method and the position of the protein
		
		r(1, :) = [0,0,0]
		
		do i=1, N-1
			call random_number(seed)
			lambda = acos(2*seed(1) - 1) - (PI/2)
			phi = 2*PI*seed(2)
			sphere = [b*cos(lambda)*cos(phi), b*cos(lambda)*sin(phi), b*sin(lambda)]
			r(i+1, :) = r(i, :) + sphere
		end do
		
		! We create an array s which will simplify the update rule in the Rouse model
		s(2:N+1, :) = r
		s(1, :) = s(2, :)
		s(N+2, :) = s(N+1, :)

		call random_number(seed)
		lambda = acos(2*seed(1) - 1) - (PI/2)
		phi = 2*PI*seed(2)
		sphere = [d0(y)*cos(lambda)*cos(phi), d0(y)*cos(lambda)*sin(phi), d0(y)*sin(lambda)]
		p = r(51, :) + sphere
		
		suml = suml + norm2(r(N,:) - r(1,:))
		
		! Now we are ready to start simulation
		
		do while(sts)
			p = p + sqrt(2*Dp*dt)*[random_normal(), random_normal(), random_normal()]
			
			do j=1, N
				rnew(j, :) = s(j+1, :) + (k/z)*(s(j, :) - 2*s(j+1, :) + s(j+2, :))*dt + &
				sqrt((2*kB*T*dt)/z)*[random_normal(), random_normal(), random_normal()]
			end do
			
			
			r = rnew
			s(2:N+1, :) = r
			s(1, :) = s(2, :)
			s(N+2, :) = s(N+1, :)
			
			if (min_distance(r, p) < db) then
				!print*, "Binding Successful"
				tot = tot + 1
				sts= .FALSE.
			else if (min_distance(r, p) > di) then
				!print*, "Binding Failed"
				tot = tot + 0
				sts=.FALSE.
			end if	

		end do
	end do

print*, tot, suml

end do

end program

function min_distance(r, p)

	implicit none
	integer, parameter :: long = selected_real_kind(15, 50)
	integer, parameter :: N=101
	real(kind=long), dimension(N, 3), intent(in) :: r
	real(kind=long), dimension(3), intent(in) :: p
	real(kind=long), dimension(N) :: distances
	real(kind=long) :: min_distance
	integer :: i
	
	do i=1, N
		distances(i) = sqrt(dot_product(p-r(i,:), p-r(i,:)))
	end do
	
	min_distance = minval(distances)
	
end function



