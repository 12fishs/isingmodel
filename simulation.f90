program ising
! TWO-DIMENSIONAL ISING MODEL: "Maybe it works, maybe it DOesn't."
!       program ising2

! INITIALIZE:
! nn*nn     -  is size 
! Temp      - is temperature in units of J/k_B 
! warm      - is # of time steps for warm-up, 
! mcs       - is total # time steps 
! ss(nn,nn) - are the spins 
! ibr       - keeps track of spin neighbors
!  Then open file(s) to write stuff on.
!  Then initial conditions for spins. 
!  Then initialize junk variables for averages.  
!  Then neighbour table to DO periodic boundary conditions.

!***********************************************************
! Defining variables

! force all variables to be declared explicitly
  implicit none

  real*8 aran,Temp,prob,energy,mag,fooe,foom,foomi
  real*8  susceptibility, foom2, capacity, fooe2, fooei, energy2, mag2
  real*8::ran3
  integer itime,nnn,ii,ix,iy,itemp
  integer datapoints
  integer :: iseed
  integer ::nn,warm,mcs
  integer,allocatable:: ss(:,:),ibr(:,:)
!
!**********************************************
! FOR TIMING
  real*8 STARTTIME,STIME0,STIMEREF, &
       STIMEK,STIMEKTOT,TMTIME0,TMTIME, &
       ETIME0,STIMEMPI,rhotime,rate,dclock
  integer*8::stime,stime1,cm,cr
  CALL system_clock(count_rate=cr)
  CALL system_clock(count_max=cm)
  

  rate = REAL(cr)
  WRITE(*,*) "system_clock rate ",rate
!  STARTTIME = DCLOCK()
!  call SYSTEM_CLOCK(stime)
!
!**********************************************
! Reading parameters from input file
  open(76,file='input.dat')
  read(76,*) nn             ! dimension of the cell
!  read(76,*) warm          ! warming up steps   
!  read(76,*) mcs            ! total number of monte carlo steps
  close(76) 
  warm = nn * nn * 100 !500 warmup steps per site
  mcs = warm + (nn * nn * 5 * 30000) ! (nn* nn* frequency of data per site * number of points per temp) 

  
!
!**********************************************
! setting dimensions of arrays
  allocate(ss(1:nn,1:nn))
  allocate(ibr(1:nn,2))
!
!**********************************************
! Opening files for writing 
!  open (77,file='configuration.dat')
  open (78,file='mag.dat')
!
!********************************************** 
!
! Initialising the arrays of spins in a 
! ferromagnetic configuration 
!
     DO  ix=1,nn 
        DO  iy=1,nn 
           ss(ix,iy) = 1 
        END DO
     END DO
!
! Iinitialising the informations about the 
! neighbours of each spin including periodic 
! boundary conditions
!
     DO ix=1,nn
        ibr(ix,1)=ix-1
        ibr(ix,2)=ix+1
        if(ix.eq.1)ibr(ix,1)=nn
        if(ix.eq.nn)ibr(ix,2)=1
     END DO
!
!*********************************************
  DO itemp=1,125           ! Loop over Temperatures

     Temp=1.5 + itemp * 0.0104

     datapoints = 0
! seeding the random number generator
     iseed=-1288833373

! initialising the total energy and total
! magnetisation to 0     
     fooe=0.0d0
     foom=0.0d0
     foom2 = 0.0d0
     fooe2 =0.0d0
!     
!*********************************************
! LOOP: "warm" monte carlo steps to warm up
!
     DO itime=1,warm
        call mcmove(ss,nn,Temp,iseed,ibr)
     END DO
!
!********************************************     
! LOOP: "mcs-warm" monte carlo steps to 
! calculate average energy and magnetization   
     DO itime=1,mcs-warm
        CALL mcmove(ss,nn,Temp,iseed,ibr)

! fooei will be the total energy at a given time, like foomi is the total moment is
      fooei=0.d0

      foomi=0.0d+00
      
      if (MOD(itime, nn * nn * 5).eq.0) then
        DO ix = 1,nn
          DO iy = 1,nn
             fooei = fooei-ss(ix,iy)*ss(ix,ibr(iy,1)) &
                  -ss(ix,iy)*ss(ibr(ix,2),iy)
             foomi = foomi + ss(ix,iy)
             
          END DO
       END DO
! the dabs forces us to make the magnetization positive, otherwise average can always be zero
       fooei = fooei / (nn * nn)
       foomi = foomi / (nn * nn)
       foom=foom+(foomi)
       fooe=fooe + fooei
! fooe and foom are now the sums of E and M over many simulation steps
       foom2 = foom2 + foomi * foomi
       fooe2 = fooe2 + fooei * fooei
! foom2 and fooe2 are the sums of **2 and M**2 over many simulation steps
       datapoints = datapoints + 1
      endif
    END DO
       
!
!*********************************************
! Normalisation of total values to number of 
! spins and number of monte-carlo steps

! now calculate <E> <M> and the averages of the squares <E**2> <M**2>

    energy = fooe/(dfloat(datapoints))
    mag = foom/(dfloat(datapoints))
    energy2 = fooe2 / (dfloat(datapoints))
    mag2= foom2 / (dfloat(datapoints))
 
! the factors of nn now turn these into energy and magnetization per spin, not total for the system
    susceptibility = (mag2-mag**2)  * (1.0 / Temp) 
    capacity = (energy2-energy**2)  * (1.0 / (Temp * Temp))

!
!*********************************************
! Generating the analytic result for the
! magnetization
!    
    if (Temp.ge.2.269) then
       foom=0.0
    else
       foom = (1.0 - (sinh (2.0/Temp) )**(-4.0) )**(1.0/8.0)
    endif
!*********************************************
    write (6,*) 'Temperature:', itemp,Temp 
    write (6,*) 'Energy/spin, Mag/spin',energy,mag
    write (6,*) 'Analytically: -2 to 0, and', foom
    write (6,*) 'Last configuration on file 77'
    write (78,'(6f16.8)') temp, energy, dabs(mag), foom, susceptibility, capacity
  END DO
!  DO ix = 1,nn
!     DO  iy = 1,nn
!        write (77,*) ix,iy,ss(ix,iy)
!     END DO
!     write(77,*)
!  END DO
!*********************************************
! Closing the files
!  
!  close(77)
  close(78)
!
!*********************************************
! releasing the arrays
!
  deallocate(ibr,ss)
!
!*********************************************
! Timing
!
!  STIME0 = DCLOCK()-STARTTIME
!  write(6,FMT=1001) STIME0
!  call SYSTEM_CLOCK(stime1)
!  write(6,FMT=1002) (stime1-stime)/rate
!
!*********************************************
! Writing formats used for data output
!1001 format (' CPU_TIME   : ',f9.2)
!1002 format (' WALL_TIME   : ',f9.2)
  
end program ising
!
!============================================================================
!
subroutine mcmove(ss,nn,Temp,iseed,ibr)

! ONE MONTE CARLO STEP by Metropolis: Flip probability 1 if Enew < Eold, 
! else prob is exp -(Enew-Eold)/T.  Simplified here since only there 
! are five cases in d=1 for external field = 0.
! FLIP WITH prob1 prob2  1.0   1.0   1.0   (Below spins called)
!             +     -     -     -     -           ss2
!            +++   +++   ++-   ++-   -+-      ss1 ss0 ss3
!             +     +     +     -     -           ss4
 
  integer nn,ix,iy,iix,iiy,ss0,ss1,ss2,ss3,ss4,de
  integer iseed,flip
  integer ix1,ix2,iy1,iy2
  integer ss(nn,nn),ibr(nn,2)
  real*8 prob1,prob2,Temp
  real*8 ran3
  prob1 = exp(-8.0/Temp)
  prob2 = exp(-4.0/Temp)
  
  DO iix=1,nn
     DO iiy=1,nn
        flip=1
        ix = 1 + nn * ran3(iseed)
        iy = 1 + nn * ran3(iseed)
        ss0=ss(ix,iy)
        ss1=ss(ibr(ix,1),iy)
        ss2=ss(ix,ibr(iy,1))
        ss3=ss(ibr(ix,2),iy)
        ss4=ss(ix,ibr(iy,2))
        de=2*ss0*(ss1+ss2+ss3+ss4)
        IF(de.eq.8.and.ran3(iseed).gt.prob1) THEN
           flip=1
        ELSE IF (de.eq.4.and.ran3(iseed).gt.prob2) THEN
           flip=1
        ELSE
           flip=-1
           
        END IF
        ss(ix,iy) = flip*ss(ix,iy)
     END DO
  END DO
  
end subroutine mcmove
!
!===============================================================
!  RANDOM # GENERATOR
function ran3(idum)
  integer*8 :: mbig,mseed,mz
  real*8::fac
  real*8::ran3
  integer*8::i,iff,ii,inext,inextp,k
  integer*8:: mj,mk
  integer*8::ma
  parameter (mbig=10000000,mseed=1618033,mz=0,fac=1.0/mbig)
  dimension ma(55)
  integer idum
  save ma
  save iff,inext,inextp
  data iff /0/
  
  if(idum.lt.0.or.iff.eq.0)then
     iff=1
     mj=abs(mseed-abs(idum))
     
     mj=mod(mj,mbig)
     ma(55)=mj
     mk=1
     DO i=1,54
        ii=mod(21*i,55)
        ma(ii)=mk
        mk=mj-mk
        if(mk.lt.mz)mk=mk+mbig
        mj=ma(ii)
     END DO
     DO k=1,4
        DO i=1,55
           ma(i)=ma(i)-ma(1+mod(i+30,55))
           if(ma(i).lt.mz)ma(i)=ma(i)+mbig
        END DO
     END DO
     inext=0
     inextp=31
     idum=1
  endif
  inext=inext+1
  if(inext.eq.56)inext=1
  inextp=inextp+1
  if(inextp.eq.56)inextp=1
  mj=ma(inext)-ma(inextp)
  if(mj.lt.mz)mj=mj+mbig
  ma(inext)=mj
  ran3=mj*fac
end function ran3
