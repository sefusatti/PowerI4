    module   input_catalog
    use      interpolation
    implicit none

    contains

      
! *****************************************************************************!
! Choose input file type. You can add here your own additional file type       !
! Particles are assigned individually by calling subroutine call assign(r,L)   !
! with double precision r(3) containing the particle coordinates and double    !
! precision L(3) containing the box sides. Also, your subroutine must count    !         ! the number of particle Nptot                                                 !
      
    subroutine input_chosen_catalog(infiletype)
        
      integer :: infiletype
      
      if(infiletype.eq.0) call input_catalog_xyz
      if(infiletype.eq.1) call input_catalog_gadget
        
    end subroutine input_chosen_catalog


      
! *****************************************************************************!
! Reading unformatted files with real, single precision coordinates x,y,z on   !
! each line. To read ascii files uncomment the read statement below            !
     
   subroutine input_catalog_xyz

     use parbox, only: Rx,Ry,Rz,infile,Nptot
     
     implicit none
     real (kind=8) :: r(3),L(3)
     real (kind=4) :: rr(3)

     L=(/Rx,Ry,Rz/)

     open(11,file=infile,status='old',form='unformatted')
     Nptot=0
     
     do while(.true.)
        read(11,end=10) rr(1),rr(2),rr(3)  ! unformatted
        r=rr
        !read(11,*,end=10) r(1),r(2),r(3) ! formatted
        !if (Nptot.eq.0) write(*,'(A,3F10.3)') ' 1st object position:', r(1),r(2),r(3)
        
        call assign(r,L)

        Nptot=Nptot+1  ! counting particles
        
     enddo

10   continue
     
   end subroutine input_catalog_xyz

!*******************************************************************!
!**** Reads Gadget unformatted comoving output files   *************!

   subroutine input_catalog_gadget 

     use parbox, only: Rx,Ry,Rz,infile,Nptot,iRSD
     
     implicit none
     real (kind=8)  :: r(3),L(3)
     integer        :: i, j, k, nc, Np

     real (kind=4) , allocatable :: rr(:,:),vv(:,:)
     integer       :: Npart(6),FlagSfr,FlagFeedback,FlagCooling,NdumFiles
     integer       :: FlagStellarAge,FlagMetals,NallHW(6),flagentrics,Nall(6)
     real          :: unused(15)
     real (kind=8) :: H,Massarr(6),afactor,redshift,Omega0,OmegaL0,hlittle,RBox
     character     :: c4*4, filename*400
     real(kind=4)  :: kpc2mpc

      kpc2mpc = 1000.d0

      Nptot=0
      k=0
      do while(.true.)
      write(c4,'(i4)') k
      nc=int(log10(float(max(k,1))))
      filename=trim(infile)//'.'//c4(4-nc:4)

      write(*,'(A,A)') ' reading: ', trim(filename)

      open(11,file=filename,status='old',err=104,form='unformatted',convert='big_endian')

      read(11) Npart,Massarr,afactor,redshift,FlagSfr,FlagFeedback,Nall,  &
               FlagCooling,NdumFiles,RBox,Omega0,OmegaL0,hlittle,         &
               FlagStellarAge,FlagMetals,NallHW,flagentrics,unused

!      write(*,*)'header read'
!      write(*,*)'Npar = '      ,  Npart(2)
!      write(*,*)'redshift z = ', redshift

      Np=Npart(2)

      Rbox = Rbox / kpc2mpc
      if(Rbox.ne.Rx) then
      write(*,'(A)') '> ERROR: Make sure your input box-size is correct'
      write(*,'(A,3F13.2)') '> Box sizes from parameters file:', Rx,Ry,Rz
      write(*,'(A,1F13.2)') '> Box size from catalog file:', Rbox
      stop
      endif 

      L=(/Rx,Ry,Rz/)

      allocate(rr(3,Np))
      allocate(vv(3,Np))
      read(11)((rr(i,j),i=1,3),j=1,Np)
      read(11)((vv(i,j),i=1,3),j=1,Np)
      close(11)
       
      H = 100.0*sqrt(Omega0*(1+redshift)**3.0+OmegaL0)/(1.0+redshift)

        do j=1,Np

           rr(1,j) = rr(1,j) / kpc2mpc   
           rr(2,j) = rr(2,j) / kpc2mpc   
           rr(3,j) = rr(3,j) / kpc2mpc   

           ! add peculiar velocity along axis iRSD
           if(iRSD.eq.1) rr(1,j) = rr(1,j) + vv(1,j)/real(H*dsqrt(1.d0+redshift))
           if(iRSD.eq.2) rr(2,j) = rr(2,j) + vv(2,j)/real(H*dsqrt(1.d0+redshift))
           if(iRSD.eq.3) rr(3,j) = rr(3,j) + vv(3,j)/real(H*dsqrt(1.d0+redshift))

           r(1) = dble(rr(1,j))
           r(2) = dble(rr(2,j))
           r(3) = dble(rr(3,j))

           call assign(r,L)

           Nptot=Nptot+1  ! counting particles        

        enddo

       deallocate(rr,vv)

       k=k + 1
       enddo


104   if(k.eq.0) then
      write(*,*) '> ERROR = file not found. Stopped',trim(filename)
      stop
      else
      write(*,'(A,I5)') ' done reading input. Number of input files',k
      endif
         
      return
      end subroutine input_catalog_gadget 

    end module
