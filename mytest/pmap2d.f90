  integer, parameter :: nx=256*2,nz=256*2

  real, dimension(nx,nz) :: u0
  character*4 fn3
  character*1 fn1
  character(100),parameter :: dir='data/r1.05neil'

  open(11,file=trim(dir)//'/icount.dat')
  read(11,*) nt
  write(*,*) 'nt=',nt
  open(12,file=trim(dir)//'/u0.dat',access='stream')
  open(13,file=trim(dir)//'/udiv.dat',access='stream')
  open(14,file=trim(dir)//'/ucrl.dat',access='stream')
  open(15,file=trim(dir)//'/ucrly.dat',access='stream')
  do ii=12,15
  write(fn1,'(i1)') ii-12
  rmax=-1.e-30
  rmin=1.e30
  do it=1,nt
  read(ii) u0
  if (ii .eq. 12) u0=u0-u0(nx,nz-1)
  rmax=max(rmax,maxval(u0))
  rmin=min(rmin,minval(u0))
  enddo
  rewind(ii)
  write(*,*) rmax,rmin
  do it=1,nt
  write(fn3,'(i4.4)')it
  read(ii) u0
  if (ii .eq. 12) u0=u0-u0(1,1)!u0(nx,nz-1)
  u0(1,1)=rmin
  u0(nx,nz)=rmax
  iscale=ii-11
  if (ii .eq. 15) iscale=2
  if (ii .eq. 12) iscale=2
  call pmap('Movie/u'//fn1//'.'//fn3//'.pgm',u0,nx,nz,iscale)
  enddo
  enddo

contains
  subroutine pmap(fn,rmap1,nx,ny,iscale)
  real rmap(nx,ny),rmap1(nx,ny),rmax,rmin
  integer*2, dimension(nx,ny) :: imap
  integer*1, dimension(nx,ny) :: imap1
  character(len=*):: fn
  integer npix,mypos

  npix=min(ny/2-1,nx/2-1,300)
  
  
  rmap=rmap1
  iscale1=iscale
  do while (iscale1 > 1)      
     rmap=sign((sqrt(abs(rmap))),rmap)
     iscale1=iscale1-1
  end do
  rmax=maxval(rmap)
  rmin=minval(rmap)
  rmap(nx/2:nx/2+1,ny/2:ny/2+1)=rmax
  write(*,*) trim(fn),rmax,rmin
  imap=225-255*(rmap-rmin)/(rmax-rmin)
  imap1=127-127*(rmap-rmin)/(rmax-rmin)
  open(10,file=fn)
  write(10,'(2hP5)')
  write(10,*)nx,ny
  write(10,*) 255
!  write(10,*) 127
  INQUIRE(UNIT=10, POS=mypos)
  close(10)
  open(10,file=fn, access='stream',position='append')
!  write(10,pos=mypos) int(imap,1)
  write(10) int(imap,1)
  close(10)
end subroutine pmap

end program
