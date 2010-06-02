c -*- Fortran -*-

c common block holding all required parameters

      logical intflg

c simulation parameters:
c ict, ictpres, lmo, akmax0, akmaxt, nlnkt:

      integer p_ict
      integer p_lmx
      real    p_emx
      real    p_kmx
      integer p_nkl

c transfer-function parameters:
c ntf, ztf, filenames:

      integer p_ntf
      real    p_ztf(10)
      character*500 p_ftf(10)

c cosmological parameters:
c omegab,omegac,omegav,omegan,h0,tcmb,yhe,annur,annunr,gsnunr:

      real    p_par(10)

c quintessence, 5dim

      integer    p_quin
      real       p_wdyn
      integer    p_dim

c recombination parameter:
c rcflag

      integer p_rcf

c reionization parameters:
c riflag, optdlss, zri, rif:

      real    p_rif
      real    p_opd
      real    p_red
      real    p_frc

c mode specifications:
c itflag, nn, an, itn, ant, irt, rat:

      integer p_tfl
      integer p_nsm
      real    p_scl(10)
      real    p_asc(10)
      real    p_das(10)
      real    p_ate(10)
      integer p_itn
      real    p_ten(10)
      integer p_irt
      real    p_rat(10)

c lensing:
c lensflag, filename for scalar C_l, filename for lensed scalar C_l,
c filejlens:

      integer p_glf
      character*500 p_fcl
      character*500 p_fgl
      character*500 p_jgl

c filename for tensor C_l:

      character*500 p_tcl

c initial-condition flag:

      integer p_inf

c filenames for j_l and u_jl:

      character*500 p_fjl
      character*500 p_ujl

c filename for FITS output:

      character*500 p_fts

      common /cmbpar/
     &     p_ict,p_lmx,p_emx,p_kmx,p_nkl,p_ntf,p_ztf,
     &     p_ftf,p_par,p_rif,p_opd,p_red,p_frc,p_tfl,p_nsm,
     &     p_scl,p_asc,p_das,p_itn,p_ten,p_ate,p_irt,p_rat,p_glf,
     &     p_fcl,p_fgl,p_jgl,p_tcl,p_inf,p_fjl,p_ujl,p_fts,p_rcf,
     &     p_quin,p_wdyn,p_dim,intflg

      save /cmbpar/

c end of parameter specification
