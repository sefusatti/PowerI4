     module grid
       integer Nx,Ny,Nz
       complex(kind=8), dimension(:,:,:), allocatable, save :: dcl
       real(kind=8), dimension(:,:,:), allocatable, save :: dtl
     end module grid
