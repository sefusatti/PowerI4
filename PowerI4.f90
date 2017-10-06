    program PowerI4

!{ Developers = Emiliano Sefusatti, Martin Crocce !}
!{ Based on previous codes by R. Scoccimarro and H. Couchman - !}
!{ Version    = beta       !} 
!{ Time-stamp = Sep - 2017 !}
!{ 2017-09-29
!{- FFT in place for interlacing
!{- Simplified indexing in assign subroutines
!{- Normalization by 1/Nptot after assignment
!{- Normalization by 1/kF3 in fcomb

!{ Subroutines in this file:
!{ input_parameters_file     
!{ print_initial_log
!{ input_parameters_file           
!{ output_real_space_density
!{ output_fourier_space_density
!{ fcomb
!{ measure_pk 

!{ Module parbox      

!{ Module grid 

!{ Module input_catalog:
!{ input_catalog_xyz     
!{ input_catalog_gadget      

!{ Module interpolate_particles:
!{ assign
!{ assign_direct
!{ assign_double_grid      
!{ assign_single_grid
!{ assign_index
!{ assign_weights_2
!{ assign_weights_3
!{ assign_weights_4


! {ciccia
      
    use input_catalog
    use parbox 
    use grid
    use, intrinsic   :: iso_c_binding ! for FFTW3       
    implicit         none
    integer(kind=8)  :: planf
      
    include 'fftw3.f'  ! FFTW3 header

!{  input parameters from external file  !}

      call input_parameters_file
      call print_initial_log 
 
!{  Allocating grids !}

      write(*,'(A)') ' allocating grid'
      if (nintp.ne.0) then
         if (interlacing) then
            ! allocate a single complex array
            allocate(dcl(Nx,Ny,Nz)) 
         else if (.not. interlacing) then
            ! allocate a complex and real arrays 
            allocate(dcl(Nx/2+1,Ny,Nz),dtl(Nx,Ny,Nz))
         endif
         rupr=(1.d0-tiny(1.d0))  !make sure there is no problem later with interp.
      elseif (nintp.eq.0) then
         ! allocate a single complex array
         allocate(dcl(Nx/2+1,Ny,Nz))
      endif

!{ Initializing grids!}

      if (.not. interlacing) dtl=0.d0 
      dcl=dcmplx(0.d0,0.d0)

!{ Reading catalogs and assigning particles !}

      write(*,'(A)') ' doing particle assignment'

      call input_chosen_catalog(infiletype)
      
      write(*,'(A,I12)') ' total number of particles = ',Nptot

!{ Box & grid parameters !}
      
      kFx=twopi/Rx             ! fundamental frequencies
      kFy=twopi/Ry
      kFz=twopi/Rz
      kF3=kFx*kFy*kFz

      kNx=kFx*dble(Nx)/2.d0    ! Nyquist frequencies 
      kNy=kFy*dble(Ny)/2.d0   
      kNz=kFz*dble(Nz)/2.d0   
       
       write(*,'(A,3F12.6)') ' fundamental frequencies = ',kFx,kFy,kFz
       write(*,'(A,3F12.6)') ' Nyquist frequencies = ',kNx,kNy,kNz

!{ Normalization !}

      if (nintp.eq.0) then
         dcl=dcl/dble(Nptot)/kF3  
      else if (interlacing) then
         dcl=dcl/dble(Nptot) 
      else
         dtl=dtl/dble(Nptot) 
      endif
      
!{ FFT !}

       if (nintp.ne.0) then 

          write(*,'(A)') ' doing FFT' 
       
          if (interlacing) then
             ! FFTW complex to complex in place 
             call dfftw_plan_dft_3d(planf,Nx,Ny,Nz,dcl,dcl,FFTW_ESTIMATE)
             call dfftw_execute(planf)
          else
             ! FFTW real to complex out of place
             call dfftw_plan_dft_r2c_3d(planf,Nx,Ny,Nz,dtl,dcl,FFTW_ESTIMATE)
             call dfftw_execute(planf)
          endif

          write(*,'(A)') ' doing fcomb'

          call fcomb

       endif
       
       call output_fourier_space_density ! only if required

       call measure_pk

       write(*,'(A)') ' all done'
       write(*,*)
       
    stop
    end program PowerI4

    
    
! *************************************************************************************
! *** Input parameters ****************************************************************

    
    subroutine input_parameters_file

    use parbox
    use grid, only: Nx,Ny,Nz
    implicit none
    character(5) :: car
       
      read(*,*) ! input file name
      read(*,*) infile
      read(*,*) ! type of input file
      read(*,*) infiletype
      read(*,*) ! output file
      read(*,*) outfile
      read(*,*) ! linear grid sizes
      read(*,*) Nx
      read(*,*) Ny
      read(*,*) Nz
      Nv=(/Nx,Ny,Nz/)
      read(*,*) ! B-spline interpolation order
      read(*,*) nintp
      if (nintp.ge.3) nintph=3
      if (nintp.lt.3) nintph=1
      read(*,*) ! interlacing
      read(*,*) car
      read(car,*) interlacing
      read(*,*) ! box sizes
      read(*,*) Rx
      read(*,*) Ry
      read(*,*) Rz
      read(*,*) ! center of the first bin in units of the fundamental frequency 
      read(*,*) rkF
      read(*,*) ! size of the bin in units of the fundamental frequency 
      read(*,*) rdk
      read(*,*) ! measure multiples
      read(*,*) iRSD
      read(*,*) ! choice output FFT density file
      read(*,*) ioutput
      ! ioutput = 1 / 2 : real-space / Fourier-space density
      ! ioutput = -1/-2 : real-space / Fourier-space density, stop after output
      read(*,*) ! output density file
      read(*,*) file_density 
       
    end subroutine input_parameters_file

    
    
! *************************************************************************************
! *** Print inputs to screen **********************************************************

    
    subroutine print_initial_log

    use parbox
    use grid, only: Nx,Ny,Nz
    implicit none

      write(*,*)
      write(*,'(A)')        ' *** PowerI4 *** '
      write(*,*)
      write(*,'(A,A)')      ' input file = ', trim(infile)
      write(*,'(A,3F10.3)') ' box sizes [Mpc/h] = ',    Rx,Ry,Rz
      write(*,'(A,3I5)')    ' FFT grid sizes = ',       Nx,Ny,Nz
      write(*,'(A,I5)')     ' interpolation order = ',  nintp
      if (Nx.ne.nint(Nz*Rx/Rz) .or. Ny.ne.nint(Nz*Ry/Rz)) then
         write(*,'(A)')     ' WARNING: not equal resolution in all directions'
      endif
      if (nintp.eq.0) then
         write(*,'(A)')     ' WARNING: direct summation!'
      endif
      if (.not. interlacing) then
         write(*,'(A)')     ' WARNING: no interlacing!'
      endif
      if(iRSD.ne.0) then
         write(*,'(A,I2)')  ' multiples evaluated along axis ', iRSD
      endif
      write(*,'(A,A)')      ' output file = ', trim(outfile)

    end subroutine print_initial_log


        
! ************************************************************************************
! *** fcomb **************************************************************************

    
    subroutine fcomb

    use parbox, only : nintp, interlacing, twopi, kF3
    use grid
    implicit none
    integer :: Nnyqx,Nnyqy,Nnyqz,icx,icy,icz,ix,iy,iz,i
    real(kind=8) :: tpiNx,piNx,tpiNy,piNy,tpiNz,piNz,cf,rk,Wkx,Wky,Wkz,cfac
    complex(kind=8) :: recx,recy,recz,xrec,yrec,zrec
    complex(kind=8) :: c1,ci,c000,c001,c010,c011,cma,cmb,cmc,cmd
    real(kind=8) :: tWkx(Nx/2+1),tWky(Nx/2+1),tWkz(Nx/2+1)

       Nnyqx=Nx/2+1
       tpiNx=twopi/dble(Nx)
       Nnyqy=Ny/2+1
       tpiNy=twopi/dble(Ny) 
       Nnyqz=Nz/2+1
       tpiNz=twopi/dble(Nz)

       tWkx(1)=1.d0
       do i=2,Nnyqx
          rk=tpiNx*dble(i-1)
          tWkx(i)=(dsin(rk/2.d0)/(rk/2.d0))**nintp
       enddo
       tWky(1)=1.d0
       do i=2,Nnyqy
          rk=tpiNy*dble(i-1)
          tWky(i)=(dsin(rk/2.d0)/(rk/2.d0))**nintp
       enddo
       tWkz(1)=1.d0
       do i=2,Nnyqz
          rk=tpiNz*dble(i-1)
          tWkz(i)=(dsin(rk/2.d0)/(rk/2.d0))**nintp
       enddo

       if (interlacing) then

          cf=1.d0/(4.d0*kF3)
          
          piNx=-tpiNx/2.d0
          recx=dcmplx(dcos(piNx),dsin(piNx))
          piNy=-tpiNy/2.d0
          recy=dcmplx(dcos(piNy),dsin(piNy))
          piNz=-tpiNz/2.d0
          recz=dcmplx(dcos(piNz),dsin(piNz))

          c1=dcmplx(1.d0,0.d0)
          ci=dcmplx(0.d0,1.d0)

          zrec=c1
          do iz=1,Nnyqz
             icz=mod(Nz-iz+1,Nz)+1
             Wkz=tWkz(iz)
             
             yrec=c1
             do iy=1,Nnyqy
                icy=mod(Ny-iy+1,Ny)+1
                Wky=tWky(iy)
                
                xrec=c1
                do ix=1,Nnyqx
                   icx=mod(Nx-ix+1,Nx)+1
                   Wkx=tWkx(ix)

                   cfac=cf/(Wkx*Wky*Wkz)
                   
                   cma=ci*xrec*yrec*zrec
                   cmb=ci*xrec*yrec*dconjg(zrec)
                   cmc=ci*xrec*dconjg(yrec)*zrec
                   cmd=ci*xrec*dconjg(yrec*zrec)
                   
                   c000=dcl(ix,iy ,iz )*(c1-cma)+dconjg(dcl(icx,icy,icz))*(c1+cma)
                   c001=dcl(ix,iy ,icz)*(c1-cmb)+dconjg(dcl(icx,icy,iz ))*(c1+cmb)
                   c010=dcl(ix,icy,iz )*(c1-cmc)+dconjg(dcl(icx,iy ,icz))*(c1+cmc)
                   c011=dcl(ix,icy,icz)*(c1-cmd)+dconjg(dcl(icx,iy ,iz ))*(c1+cmd)
                 
                   dcl(ix,iy ,iz )=c000*cfac
                   dcl(ix,iy ,icz)=c001*cfac
                   dcl(ix,icy,iz )=c010*cfac
                   dcl(ix,icy,icz)=c011*cfac
                   dcl(icx,iy ,iz )=dconjg(dcl(ix,icy,icz))
                   dcl(icx,iy ,icz)=dconjg(dcl(ix,icy,iz ))
                   dcl(icx,icy,iz )=dconjg(dcl(ix,iy ,icz))
                   dcl(icx,icy,icz)=dconjg(dcl(ix,iy ,iz ))
                 
                   xrec=xrec*recx
                enddo
                yrec=yrec*recy
             enddo
             zrec=zrec*recz
          enddo

     else ! no interlacing, only correcting for the window function and normalizing

        cf=1.d0/kF3

        do iz=1,Nnyqz
           icz=mod(Nz-iz+1,Nz)+1
           Wkz=tWkz(iz)
              
           do iy=1,Nnyqy
              icy=mod(Ny-iy+1,Ny)+1
              Wky=tWky(iy)
              
              do ix=1,Nnyqx
                 Wkx=tWkx(ix)
                 
                 cfac=cf/(Wkx*Wky*Wkz)                 
                 dcl(ix,iy ,iz )=dcl(ix,iy ,iz )*cfac
                 if(iz.ne.icz) dcl(ix,iy,icz)=dcl(ix,iy,icz)*cfac
                 if(iy.ne.icy) dcl(ix,icy,iz)=dcl(ix,icy,iz)*cfac
                 if(iz.ne.icz .and. iy.ne.icy) dcl(ix,icy,icz)=dcl(ix,icy,icz)*cfac

              enddo
           enddo
        enddo        

     endif
         
    end subroutine fcomb


    
! ***********************************************************************************
! *** Output densities  *************************************************************

  
    
    subroutine output_fourier_space_density
      
      use parbox, only : ioutput,file_density,Nptot,kFx,kFy,kFz,kNx
      use grid, only   : Nx,Ny,Nz,dcl
      implicit none
      integer :: ix,iy,iz
      
      if (abs(ioutput).eq.1) then
         open(unit=4,file=file_density,status='unknown',form='unformatted')
         write(4) Nx,Ny,Nz,dble(Nptot),kFx,kFy,kFz,kNx
         write(4) (((dcl(ix,iy,iz),ix=1,Nx/2+1),iy=1,Ny),iz=1,Nz)
         close(4)
         if (ioutput.lt.0) then
            write(*,'(A)') ' all done '
            write(*,*)
            stop
         endif
      endif

      
    end subroutine output_fourier_space_density


        
! ***********************************************************************************
! *** Measure and output P(k)  ******************************************************    

    
    subroutine measure_pk

      use parbox
      use grid, only: Nx,Ny,Nz,dcl
      implicit none
      integer :: ix, iy, iz, ikx, iky, ikz, icx, icy, icz, Nbins, i, imk
      real (kind=8), dimension(:), allocatable :: avgk, avgP, avgP2, avgP4, co
      real (kind=8) :: k,kFmin,rk,pSN,kNmax,dk
      real (kind=8) :: Le2, Le4, sp, sit1, cot1, rkx, rky, rkz, pk, coga, cp, cc
      real (kind=8) :: thetaobs, cthetaobs, sthetaobs, phiobs
    complex (kind=8) :: ct  

      write(*,'(A)') ' measuring power spectrum'
       
      kFmin=min(kFx,min(kFy,kFz))
      kNmax=max(kNx,max(kNy,kNz)) 

      ! bin size is determined based on the smallest fundamental frequency
      ! while the number of bins depends on the largest Nyquist frequency 

      dk=rdk*kFmin   ! size of the wavenumbers bin
      Nbins=nint(kNmax/dk)
      write(*,'(A,I5)') ' number of bins = ', Nbins
       
      pSN=1.d0/(kF3*dble(Nptot))     ! shot-noise contribution for Poisson distribution
       
      allocate (avgk(Nbins),avgP(Nbins),co(Nbins))
      avgk=0.d0
      avgP=0.d0
      co=0.d0

      thetaobs=0.d0          
      sthetaobs=dsin(thetaobs)
      cthetaobs=dcos(thetaobs)
      phiobs=0.d0 
      if (iRSD.eq.0) then   ! no multiples
         thetaobs=0.d0          
         phiobs=0.d0 
      else                  ! multiples
         allocate (avgP2(Nbins),avgP4(Nbins))  
         avgP2=0.d0
         avgP4=0.d0
         if (iRSD.eq.1) then
             thetaobs= 0.5d0*pi     !0.837*pi
             phiobs=0.d0            !0.595*pi
         elseif (iRSD.eq.2) then
             thetaobs=0.5d0*pi      !0.527*pi
             phiobs=0.5d0*pi        !1.53*pi
         elseif (iRSD.eq.3) then
             thetaobs=0.d0          !0.235*pi
             phiobs=0.d0            !0.783*pi
         endif
          cthetaobs=cos(thetaobs)
          sthetaobs=sin(thetaobs)
      endif

      do iz=1,Nz
         ikz=mod(iz+Nz/2-2,Nz)-Nz/2+1
         icz=mod(Nz+1-iz,Nz)+1
         rkz=kFz*float(ikz)

      do iy=1,Ny
         iky=mod(iy+Ny/2-2,Ny)-Ny/2+1
         icy=mod(Ny+1-iy,Ny)+1
         rky=kFy*float(iky)

      do ix=1,Nx/2+1
         ikx=mod(ix+Nx/2-2,Nx)-Nx/2+1
         icx=mod(Nx+1-ix,Nx)+1
         rkx=kFx*float(ikx)

         if (ix.ne.1 .or. (ix.eq.1 .and. iky.gt.0) .or. (ix.eq.1 .and. iky.eq.0 .and. ikz.gt.0)) then
   
         rk=dsqrt(rkx**2+rky**2+rkz**2)
         imk=nint(rk/dk-rkF/rdk-tiny(kFmin))+1
                
           if (imk.le.Nbins .and. imk.gt.0) then

             co(imk)=co(imk)+1.d0
             avgk(imk)=avgk(imk)+rk
             if (ix.le.Nx/2+1) then
             ct=dcl(ix,iy,iz)
             else
             ct=dcl(icx,icy,icz)
             endif
             pk=(cdabs(ct))**2
             avgP(imk)=avgP(imk)+pk
                   
             if (iRSD.ne.0) then
             cot1=rkz/rk
             sit1=dsqrt(1.d0 - cot1*cot1)
             if (sit1.gt.0.d0) then
             cp=rkx/(rk*sit1)
             sp=rky/(rk*sit1)
             cc=dsin(phiobs)*sp+dcos(phiobs)*cp
             else
             cc=0.d0
             endif
             coga = cthetaobs*cot1 + sthetaobs*sit1*cc
             Le2 = -0.5d0 + 1.5d0 * coga**2
             Le4 = 0.375d0 - 3.75d0 * coga**2 + 4.375d0 * coga**4
             avgP2(imk) = avgP2(imk) + pk*5.d0*Le2
             avgP4(imk) = avgP4(imk) + pk*9.d0*Le4
             endif

           endif
         endif
      enddo
      enddo
      enddo

      write(*,'(A)') ' writing output file'

      open(4,file=outfile,status='unknown',form='formatted')
      write(4,*) int(Nptot), pSN
      do i=1,Nbins
         if (co(i).gt.0.) then
            k = (i-1)*dk + rkF*kFmin
            avgk(i)=avgk(i) /co(i)
            avgP(i)=avgP(i) * kF3 /co(i)
            if (iRSD.ne.0) then
            avgP2(i) = avgP2(i) * kF3 /co(i)
            avgP4(i) = avgP4(i) * kF3 /co(i)
            write(4,'(5E18.8,I16)') k, avgk(i), avgP(i), avgP2(i), avgP4(i), int(co(i))
            else
            write(4,'(3E18.8,I16)') k, avgk(i), avgP(i), int(co(i))
            endif
         end if
      enddo
      close(4)

    end subroutine measure_pk    






