c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mechanic  --  initialize molecular mechanics  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mechanic" sets up needed parameters for the potential energy
c     calculation and reads in many of the user selectable options
c
c
      subroutine mechanic
      use inform
      use iounit
      use limits
      use potent
      use vdwpot
      implicit none
      real values(2)
      real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10
      integer nt
      real st1, st2, st3, st4, st5, st6, st7, st8, st9, st10
      data nt /0/, st1 /0/, st2 /0/, st3 /0/, st4 /0/, st5 /0/
      data st6 /0/, st7 /0/, st8 /0/, st9 /0/, st10 /0/
      save nt, st1, st2, st3, st4, st5, st6, st7, st8, st9, st10

c
c
c     set the bonded connectivity lists and active atoms
c
      call attach
      call active
c
c     find bonds, angles, torsions, bitorsions and small rings
c
      call bonds
      call angles
      call torsions
      call bitors
      call rings

c
c     get the base force field from parameter file and keyfile
c
c     call dtime(values, t1)
      call field
c     call dtime(values, t2)
c
c     find unit cell type, lattice parameters and cutoff values
c
      call unitcell
      call lattice
      call polymer
      call cutoffs
c
c     setup needed for potential energy smoothing methods
c
      call flatten
c
c     assign atom types, classes and other atomic information
c
      call katom
c     call dtime(values, t3)
c
c     assign atoms to molecules and set the atom groups
c
      call molecule
      call cluster
c
c     find any pisystem atoms, bonds and torsional angles
c
      call orbital
c
c     assign bond, angle and cross term potential parameters
c
      call kbond
      call kangle
      call kstrbnd
      call kurey
      call kangang
c
c     assign out-of-plane deformation potential parameters
c
      call kopbend
      call kopdist
      call kimprop
      call kimptor
c
c     assign torsion and torsion cross term potential parameters
c
      call ktors
      call kpitors
      call kstrtor
      call kangtor
      call ktortor
c
c     assign electrostatic interaction potential parameters
c
c     call dtime(values, t4)
      call kcharge
      call kdipole
      call kmpole
      call kpolar
      call kchgtrn
      call kchgflx
c     call dtime(values, t5)
c
c     assign van der Waals, repulsion and dispersion parameters
c
      call kvdw
c     call dtime(values, t6)
      call krepel
c     call dtime(values, t7)
      call kdisp
c     call dtime(values, t8)
c
c     assign solvation, metal, pisystem and restraint parameters
c
      call ksolv
      call kmetal
      call korbit
      call kgeom
      call kextra
c
c     assign electrostatic and dispersion Ewald sum parameters
c
      call kewald
c
c     set any holonomic interatomic distance constraints
c
      call shakeup
c     call dtime(values, t9)
c
c     set hybrid parameter values for free energy perturbation
c
      call mutate
c     call dtime(values, t10)
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
c     write (*,*) 'mech', st2/nt, st3/nt, st4/nt, st5/nt, st6/nt,
c    &                    st7/nt, st8/nt, st9/nt, st10/nt
c
c     quit if essential parameter information is missing
c
      if (abort) then
         write (iout,10)
   10    format (/,' MECHANIC  --  Some Required Potential Energy',
     &              ' Parameters are Undefined')
         call fatal
      end if
      return
      end
