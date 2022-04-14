      function foo() bind ( C, name="foo" )
      use atoms
      use boxes
      use files
      use inform
      use iounit
      use output
      implicit none
      real*8 foo
      integer i,j,ixyz
      integer frame,nold
      integer size
      integer freeunit
      integer trimtext
      integer, allocatable :: list(:)
      real*8 energy
      real*8, allocatable :: derivs(:,:)
      logical dosystem,doparam
      logical doenergy,doatom
      logical dolarge,dodetail
      logical domoment,dovirial
      logical doconect,dosave
      logical exist
      logical, allocatable :: active(:)
      character*1 letter
      character*240 record
      character*240 string
      character*240 xyzfile

      real values(2)
      real t1, t2, t3, t4, t5, t6, t7, t8, t9
      integer nt
      real st1, st2, st3, st4, st5, st6, st7, st8, st9, st10
      data nt /0/, st1 /0/, st2 /0/, st3 /0/, st4 /0/, st5 /0/
      data st6 /0/, st7 /0/, st8 /0/, st9 /0/, st10 /0/
      save nt, st1, st2, st3, st4, st5, st6, st7, st8, st9, st10
c
c
c     set up the structure and mechanics calculation
c
      call etime(values, t1)
      call initial
      call etime(values, t2)
      call libgetxyz
      call etime(values, t3)
      call mechanic
      call etime(values, t4)
c
c     get the desired types of analysis to be performed
c
      string = "E"
      exist = .true.
c
c     set option control flags based desired analysis types
c
      dosystem = .false.
      doparam = .false.
      doenergy = .true.
      doatom = .false.
      dolarge = .false.
      dodetail = .false.
      domoment = .false.
      dovirial = .false.
      doconect = .false.
      call upcase (string)
c
c     set option control flag to save forces or induced dipoles
c
      dosave = .false.
      call optinit
      call etime(values, t5)
      if (frcsave .or. uindsave)  dosave = .true.
c
c     perform dynamic allocation of some local arrays
c
      size = 40
      allocate (list(size))
      allocate (active(n))
      debug = .false.
c
c     reopen the coordinates file and read the first structure
c
      frame = 0
      ixyz = freeunit ()
      xyzfile = filename
      call suffix (xyzfile,'xyz','old')
      call etime(values, t6)
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
      call etime(values, t7)
c
c     decide whether to perform analysis of individual frames
c
      abort = .false.
c
c     perform analysis for each successive coordinate structure
c
      frame = frame + 1
c
c     make the call to compute the potential energy
c
      call analysis (energy)
      call etime(values, t8)
c
c     attempt to read next structure from the coordinate file
c
c
c     perform deallocation of some local arrays
c
      deallocate (list)
      deallocate (active)
c
c     perform any final tasks before program exit
c
      close (unit=ixyz)
      call final
      foo = energy
      call etime(values, t9)
      st2 = st2 + t2 - t1
      st3 = st3 + t3 - t2
      st4 = st4 + t4 - t3
      st5 = st5 + t5 - t4
      st6 = st6 + t6 - t5
      st7 = st7 + t7 - t6
      st8 = st8 + t8 - t7
      st9 = st9 + t9 - t8
      nt = nt + 1
      if (mod(nt,1000) == 0) then
      write (*,*) 'foo', st2/nt, st3/nt, st4/nt, st5/nt, st6/nt, st7/nt,
     &                   st8/nt, st9/nt
      end if
      end function

