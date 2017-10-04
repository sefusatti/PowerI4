    module   interpolation
  
    implicit none
    
    contains

!{ Subroutines:
!{ assign
!{ assign_direct
!{ assign_double_grid      
!{ assign_single_grid
!{ assign_index
!{ assign_weights_2
!{ assign_weights_3
!{ assign_weights_4


      
!*******************************************************************
!*** Assignment choice *********************************************

      
    subroutine assign(r,L)

    use parbox, only: interlacing, nintp, rupr
    implicit none
    real (kind=8) :: r(3), L(3)
    
    r(:)=mod(r(:)+L(:),L(:)) ! make sure particle are in the box
    r(:)=min(rupr,r(:)/L(:)) ! normalize to interval [0,1)

    if (nintp.eq.0) then
       call assign_direct(r)  ! direct summation
    else if (interlacing) then
       call assign_double_grid(r)  ! double grid assignment 
    else if (.not. interlacing) then
       call assign_single_grid(r)   ! single grid assignment
    endif

    end subroutine assign


    
!*******************************************************************
!*** Direct assignment *********************************************

    
    subroutine assign_direct(r)

    use grid, only : Nx,Ny,Nz,dcl
    use parbox, only : twopi
    implicit none
    integer :: ix,iy,iz,icy,icz,ikx,iky,ikz,Nnyqx,Nnyqy,Nnyqz
    real(kind=8) :: rkx,rky,rkz,product,sp,cp,rtpi(3),r(3)
      
      Nnyqx=Nx/2+1
      Nnyqy=Ny/2+1
      Nnyqz=Nz/2+1
      rtpi=r*twopi
       
      do iz=1,Nnyqz
         ikz=mod(iz+Nz/2-2,Nz)-Nz/2+1
         rkz=dble(ikz)*rtpi(3)
         icz=mod(Nz-iz+1,Nz)+1
          
         do iy=1,Nnyqy
            iky=mod(iy+Ny/2-2,Ny)-Ny/2+1
            rky=dble(iky)*rtpi(2)
            icy=mod(Ny-iy+1,Ny)+1

            do ix=1,Nnyqx
               ikx=mod(ix+Nx/2-2,Nx)-Nx/2+1
               rkx=dble(ikx)*rtpi(1)
               
               product=rkx+rky+rkz
               cp=dcos(product)
               sp=dsin(product)
               dcl(ix,iy ,iz ) = dcl(ix,iy ,iz ) + dcmplx(cp,sp)

               if (iz.ne.Nnyqz .and. iz.ne.1) then                   
                  product=rkx+rky-rkz
                  cp=dcos(product)
                  sp=dsin(product)
                  dcl(ix,iy ,icz) = dcl(ix,iy ,icz) + dcmplx(cp,sp)
               endif
                   
               if (iy.ne.Nnyqy .and. iy.ne.1) then                   
                  product=rkx-rky+rkz
                  cp=dcos(product)
                  sp=dsin(product)
                  dcl(ix,icy,iz ) = dcl(ix,icy,iz ) + dcmplx(cp,sp)
               endif
                   
               if (iz.ne.Nnyqz .and. iy.ne.Nnyqy  .and. iz.ne.1 .and. iy.ne.1) then                   
                  product=rkx-rky-rkz
                  cp=dcos(product)
                  sp=dsin(product)
                  dcl(ix,icy,icz) = dcl(ix,icy,icz) + dcmplx(cp,sp)      
               endif
                
            enddo
         enddo
      enddo
         
    end subroutine assign_direct


    
! *************************************************************************************
! *** Assign density on a double grid for interlacing *********************************

    
    subroutine assign_double_grid(r)

    use grid, only : dcl
    use parbox, only : Nv, nintp
    implicit none
    integer :: ix,iy,iz,iv(3,4),jv(3,4)
    real(kind=8) :: r(3),w1(3,4),w2(3,4),vr(3),vt(3)
    complex(kind=8) :: cpx
              
     vr(:)=dble(Nv(:))*r(:)+1.d0
     vt(:)=vr(:)+0.5d0
     
     if (nintp.eq.4) then
        
       call assign_index(int(vr),iv)
       call assign_weights_4(vr,w1)
       call assign_index(int(vt),jv)
       call assign_weights_4(vt,w2)

     elseif (nintp.eq.3) then

       call assign_index(nint(vr),iv)
       call assign_weights_3(vr,w1)
       call assign_index(nint(vt),jv)
       call assign_weights_3(vt,w2)
        
     elseif (nintp.eq.2) then

       call assign_index(int(vr),iv)
       call assign_weights_2(vr,w1)
       call assign_index(int(vt),jv)
       call assign_weights_2(vt,w2)
        
     elseif (nintp.eq.1) then

       w1=1.d0
       w2=1.d0
       iv(:,1)=mod(nint(vr(:))-1+Nv(:),Nv(:))+1
       jv(:,1)=mod(nint(vt(:))-1+Nv(:),Nv(:))+1

     endif

     do ix=1,nintp
      do iy=1,nintp
       do iz=1,nintp

          cpx=dcmplx(w1(1,ix)*w1(2,iy)*w1(3,iz),0.d0)
          dcl(iv(1,ix),iv(2,iy),iv(3,iz))=dcl(iv(1,ix),iv(2,iy),iv(3,iz))+cpx

          cpx=dcmplx(0.d0,w2(1,ix)*w2(2,iy)*w2(3,iz))
          dcl(jv(1,ix),jv(2,iy),jv(3,iz))=dcl(jv(1,ix),jv(2,iy),jv(3,iz))+cpx

       enddo
      enddo
     enddo
   
    end subroutine assign_double_grid


    
! ************************************************************************************
! *** Assign density on a single grid *************************************************

    
    subroutine assign_single_grid(r)

    use grid, only : dtl
    use parbox, only : Nv, nintp
    implicit none
    integer :: ix,iy,iz,ivr(3),iv(3,4)
    real(kind=8) :: r(3),w(3,4),vr(3)

     w=0.d0
     vr=dble(Nv)*r+1.d0

     if (nintp.eq.4) then ! PCS

       ivr=int(vr)
       call assign_index(ivr,iv)
       call assign_weights_4(vr,w)

     elseif (nintp.eq.3) then ! TSC

       ivr=nint(vr)
       call assign_index(ivr,iv)
       call assign_weights_3(vr,w)
        
     elseif (nintp.eq.2) then ! CIC

       ivr=int(vr)
       call assign_index(ivr,iv)
       call assign_weights_2(vr,w)
        
     elseif (nintp.eq.1) then ! NGP
          
       w=1.d0
       ivr=nint(vr)
       iv(:,1)=mod(ivr(:)-1+Nv(:),Nv(:))+1

     endif

     do ix=1,nintp
      do iy=1,nintp
       do iz=1,nintp
       dtl(iv(1,ix),iv(2,iy),iv(3,iz))=dtl(iv(1,ix),iv(2,iy),iv(3,iz))+w(1,ix)*w(2,iy)*w(3,iz)
       enddo
      enddo
     enddo
       
    end subroutine assign_single_grid

   

! ************************************************************************************
! *** Indices and weights ************************************************************

    
    subroutine assign_index(iv,idx)
    use parbox, only: Nv,nintp,nintph
    implicit none
    integer :: j,iv(3),idx(3,4)
     do j=1,nintp
        idx(:,j)=mod(iv(:)-nintph+j+Nv(:),Nv(:))+1
     enddo
    end subroutine assign_index

    
    subroutine assign_weights_4(v,w)   
    implicit none
    integer :: i
    real(kind=8) :: h,h2,v(3),w(3,4)
     do i=1,3
        h=v(i)-int(v(i))
        h2=h*h
        w(i,1)=(1.d0-h)**3/6.d0
        w(i,2)=4.d0/6.d0+(0.5d0*h-1.d0)*h2
        w(i,4)=h2*h/6.d0
        w(i,3)=1.d0-w(i,1)-w(i,2)-w(i,4)
     enddo
    end subroutine assign_weights_4

    
    subroutine assign_weights_3(v,w)   
    implicit none
    integer :: i
    real(kind=8) :: h,h2,v(3),w(3,4)
     do i=1,3
        h=v(i)-nint(v(i))
        h2=h*h
        w(i,1)=0.5d0*(0.5d0-h)**2
        w(i,2)=0.75d0-h2
        w(i,3)=1.d0-w(i,1)-w(i,2)
     enddo
    end subroutine assign_weights_3

    
    subroutine assign_weights_2(v,w)   
    implicit none
    integer :: i
    real(kind=8) :: h,v(3),w(3,4)
     do i=1,3
        h=v(i)-int(v(i))
        w(i,1)=1.d0-h
        w(i,2)=h
     enddo
    end subroutine assign_weights_2

   
! ************************************************************************************

  end module interpolation
