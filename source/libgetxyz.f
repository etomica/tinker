c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine getxyz  --  get Cartesian coordinate structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "getxyz" asks for a Cartesian coordinate file name,
c     then reads in the coordinates file
c
c
      subroutine libgetxyz
      use inform
      use iounit
      use output
      implicit none
      integer ixyz
      integer freeunit
      logical exist
      character*240 xyzfile

      real values(2)
      real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10
      integer nt
      real st1, st2, st3, st4, st5, st6, st7, st8, st9, st10
      data nt /0/, st1 /0/, st2 /0/, st3 /0/, st4 /0/, st5 /0/
      data st6 /0/, st7 /0/, st8 /0/, st9 /0/, st10 /0/
      save nt, st1, st2, st3, st4, st5, st6, st7, st8, st9, st10
c
c
c     try to get a filename from the command line arguments
c
c     call :time(values, t1)
      xyzfile = "molecules.xyz"
      call basefile (xyzfile)
      call suffix (xyzfile,'xyz','old')
      inquire (file=xyzfile,exist=exist)
c
c     ask for the user specified input structure filename
c
      if (.not. exist)  call fatal
c     call dtime(values, t2)
c
c     first open and then read the Cartesian coordinates file
c
      coordtype = 'CARTESIAN'
      ixyz = freeunit ()
      open (unit=ixyz,file=xyzfile,status='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
      close (unit=ixyz)
c     call dtime(values, t3)
c
c     quit if the Cartesian coordinates file contains no atoms
c
      if (abort) then
         write (iout,30)
   30    format (/,' GETXYZ  --  Cartesian Coordinates File',
     &              ' does not Contain Any Atoms')
         call fatal
      end if
c     call dtime(values, t4)
      st2 = st2 + t2
      st3 = st3 + t3
      st4 = st4 + t4
      st5 = st5 + t5
      st6 = st6 + t6
      st7 = st7 + t7
      st8 = st8 + t8
      st9 = st9 + t9
      st10 = st10 + t10
      nt = nt + 1
c     write (*,*) 'libgetxyz', st2/nt, st3/nt, st4/nt, st5/nt
      return
      end
