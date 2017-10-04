 module parbox
   real(kind=8)    :: Rx,Ry,Rz,kFx,kFy,kFz,kF3,kNx,kNy,kNz,rkF,rdk,rupr 
   integer         :: Nv(3), iRSD, nintp, nintph,ioutput, infiletype
   integer(kind=8) :: Nptot
   character*255   :: infile, outfile, file_density
   logical         :: interlacing
   real(kind=8), parameter ::  pi = 3.1415926535897932d0, twopi = 2.d0*pi      
   save
 end module parbox
