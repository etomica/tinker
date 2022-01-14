c
c
c     #########################################################
c     ##  COPYRIGHT (C) 2021 by Moses Chung & Jay W. Ponder  ##
c     ##                   All Rights Reserved               ##
c     #########################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module kexnd  --  exchange induction forcefield parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     peisiz     Pauli repulsion size value for each atom class
c     peidmp     alpha Pauli repulsion parameter for each atom class
c     peiele     number of valence electrons for each atom class
c
c
      module kexnd
      implicit none
      real*8, allocatable :: peisiz(:)
      real*8, allocatable :: peidmp(:)
      real*8, allocatable :: peiele(:)
      save
      end
