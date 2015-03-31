c     for rosetta only  bloss62 is used 
      implicit real*8(a-h,o-z)
      integer maxa,nwt,newwt,maxp,natom,numdel,
     &        ntyp,maxit,maxres,frgwind,inb,chiral
      logical convg,gly_flag
      real*8 cbox,cx,cy,cz,vdwradi(5),pi,pi3,theta,rcut,sgcut,rd
      parameter(maxa=10000,cbox=18.856,maxit=100,inb=25,
     &          ntyp=20,natyp=8,maxres=10000,max3=3*maxa)
      integer resnum(maxa),restyp(maxres),pol(20),ikpp(10,100),
     &        ihflg(maxa),rnum,ncel,iused(ntyp),imut(maxres),
     &        itgt(maxres),iatyp(maxa),iu(natyp,natyp),rnumme,
     &        map(ntyp),ind1(maxa),ind2(maxa),ibk(maxres,6),
     &        restyp1(maxres),isec(maxres,6),isecq(maxres)
     &     ,typhold(25,9,maxres),imax1(maxres),poshld(25,maxres)
     &     ,imax2(maxres),imax3(maxres)
     &     ,l,l2,l3,n,i,j,k,id,id2,id3,id4,id5,id6,numid
     &     ,nfil,ips,ips2,ipstemp(10000),o
      real*8 xp_n(maxa),yp_n(maxa),zp_n(maxa),
     &       xp_ca(maxa),yp_ca(maxa),zp_ca(maxa),
     &       xp_cb(maxa),yp_cb(maxa),zp_cb(maxa),
     &       xp_c(maxa),yp_c(maxa),zp_c(maxa),
     &       xp_o(maxa),yp_o(maxa),zp_o(maxa),
     &       xp(20000),yp(20000),zp(20000),
     &       ddg(maxres),resdep(maxres),  
     &        hd(20),asa(20),qq(maxa),rco(maxres),ele(maxa),aa,
     &        profile(maxres),ppp(maxres),psi0(maxres),phi0(maxres),
     &        prf1(maxres),ang1(3,maxres),hold(25,9,maxres),
     &   holdsec(25,3,10,maxres),tt(3),bb(3,3),fa(20,maxres)
     &   ,fa2(20,maxres),fa3(20,maxres)
     &   ,fa2_60(20,2000),fa2_65(20,2000),fa2_80(20,2000)
     &   ,xx(150),yy(150),xp0(20000),yp0(20000),zp0(20000),
     &   mydep(maxres),
     &   idt1_gcp,idt1_nongcp,
     &   idt2_gcp,idt2_nongcp,
     &   idt3_gcp,idt3_nongcp,
     &   idt1,idt2,idt3,
     &   idt2_60,idt2_65,idt2_80,
     &   rnumme_gcp,rnumme_nongcp,alp
      real*8 atmdep(maxa),hvnum(ntyp)  
      real*8 sx,sy,sz,tx,ty,tz,xm(3,4),dep(maxa,maxit),avedep,
     &    entropy(ntyp),volume(ntyp),aenv,para(ntyp),
     &       para1(ntyp),alf(natyp,natyp),bkbn(natyp,natyp),
     &       chg(ntyp,20),rmsd(9000),cntn(20),xh(4),yh(4),zh(4),
     &       ttlcnt(14),xmol(20,20),xmaxcnt(20),avecnt(20),
     &     hphi(20),hpsi(20),psec(maxres,3),
     &     bls(ntyp,ntyp),bls0(ntyp,ntyp),cov
     &     ,xpmean,ypmean,zpmean,xptmp,yptmp,zptmp,temp,tmp,tmp1
     &     radius
      real bond, angle1,angle2
      real x1(3),x2(3),x3(3),x4(3),x5(3),x6(3),cb(3),dif(4,maxres)
      character resname*4,atmname*4,ctmp*2,resn(ntyp)*3,rname(maxa)*4,
     &    aname(maxa)*4,line1*80,cct*2,chain*1,ares1*5,
     & ares_n*5, ares_ca*5, ares_cb*5, ares_c*5, ares_o*5,
     & ares_n_before*5, ares_n_before2*5, ares_n_before3*5,
     & ares_ca_before*5, ares_ca_before2*5,
     & ares_cb_before*5, ares_cb_before2*5, 
     & ares_c_before*5, ares_c_before2*5,
     & ares_o_before*5, ares_o_before2*5,
     &    ares0*5,   aresbefore*5,
     &    sshold(200,10,maxres)*1,ssnt(maxres)*1
     &    ,ssholdsec(200,10,maxres)*1
      character resnm*3,atnm*3,cind(ntyp,20)*3,hlc(25,maxres)*10,
     &   myseq(maxres)*3,hlcsec(25,maxres)*10,resmap(20)*1

        real*8 disp,disp0(25,maxres),disp0sec(25,maxres)
     &  ,disp1(25,maxres)

       data qq/maxa*0.0/,ele/maxa*0.0/,aa/0.0/
       data map/1,12,11,20,9,10,14,16,15,13,17,18,
     &          19,2,4,5,6,3,8,7/ 

       data asa/360.,425.,429.,402.,399.,375.,470.,447.,347.,
     &  323,370.,352.,408.,379.,405.,378.,432.,468.,445.,364./

       data hd/4.19,4.80,6.03,5.87,5.89,4.51,7.14,4.90,2.03,1.01,
     &     2.53,1.59,2.48,1.40,1.78,1.09,2.95,1.78,2.01,3.51/



      data hvnum/2,4,7,4,4,3,10,8,1,0,
     &           3,2,5,4,5,4,6,7,5,3/               

      data pol/1,0,0,0,0,0,0,1,0,0,1,1,1,1,1,1,1,1,1,0/

      data volume/45.,104.,130.,101.,101.,75.,168.,133.,26.,0.,
     &     56.,30.0,86.0,64.,77.,53.,96.0,129.,106.,59./



      data entropy/-0.55,-1.61,-0.58,-0.89,-0.78,-0.51,-0.97,-0.98,
     & 0.0,0.0,-1.63,-1.71,-2.11,-1.57,-1.80,-1.25,-0.96,-2.03,
     & -1.94,0.0/            

      data vdwradi,pi/2.0,1.85,1.70,2.0,0.,3.1415926/

      data resn/'CYS','MET','PHE','ILE','LEU','VAL','TRP','TYR',
     &          'ALA','GLY','THR','SER','GLN','ASN','GLU','ASP',
     &          'HIS','ARG','LYS','PRO'/
   

      data resmap/'C','M','F','I','L','V','W','Y','A','G',
     &            'T','S','Q','N','E','D','H','R','K','P'/

c     
      character whole*80, filename*80,whole1(maxa)*80,fil1*80,fil2*80,
     &           tempafil(10000)*80,afil(10000)*80,fbase*80,file0*80,
     &           firsts, targets

         call getarg(1,fbase)
         call getarg(2,file0)
!           read(*, '(a80)')  fbase
         ips=index(fbase,' ')-1
         fbase=trim(fbase)
         file0=trim(file0)
         ips2=index(file0,' ')-1

        open(11,file=fbase(1:ips)//'/lib-bin/testpdb.lst2',status='old')
!        print*, fbase
!        print*, file0

        i=1
        afil(1)=file0(1:ips2)
        i=i+1

c        open(unit=20,file='pdbfile.list',status='old')
 2001   read(11,*,end=2000) tempafil(i)
        tempafil(i)=trim(tempafil(i))
        ipstemp(i)=index(tempafil(i),' ')-1
        afil(i)=tempafil(i)(1:ipstemp(i))
        i=i+1
        goto 2001

 2000   continue
        nfil=i-1
        close(11)
c        close(20)

!      do ii=1,nfil
!        print*, afil(1)
!        print*, afil(1)(1:ips2), afil(ii)(1:ipstemp(ii)),'end'
!      enddo


      convg = .FALSE.
      bond=1.54
      angle1=109.53
      angle2=123.00
      chiral=1


       frgwind=5
       dd=1./2.8
c       cov=15.0
       cov=0.0
       alp=0.48

          do i=1,10000
            do j=1,20
               fa(j,i)=0.0
               fa2(j,i)=0.0
               fa3(j,i)=0.0
            enddo
          enddo

c        goto 3455

          do i=1,inb
               do j=1,maxres
               disp0(i,j)=1.d8
               enddo
          enddo

          do i=1,20000
              xp0(i)=0.0
              yp0(i)=0.0
              zp0(i)=0.0
          enddo
          do i=1,10000
              restyp1(i) = 0
              mydep(i) = 0.0
          enddo

          do i=1,2000
            do j=1,20
               fa2_60(j,i)=0.0
               fa2_65(j,i)=0.0
               fa2_80(j,i)=0.0
            enddo
          enddo


        jjj = 0
      do ii=1,nfil
!        if(ii .eq. 2) goto 500
!            firsts=afil(1)(1:ips2)
!            targets=afil(ii)(1:ipstemp(ii))
!        print*, firsts, targets,'end'
!         if(afil(1) .eq. '4ubp.pdb' .and. afil(ii) .eq. '1ay7b')then
!             print*, afil(ii),'Yes'
!         else 
!             print*, afil(ii),'No'
!         endif
        
         if ( (ii .ne. 1) .and. 
     &      (
     &( afil(1)(1:4) .eq. '1a19' .and. afil(ii)(1:5) .eq. '1ay7b').or.  
     &( afil(1)(1:4) .eq. '1aiu' .and. afil(ii)(1:5) .eq. '1gh2a').or.  
     &( afil(1)(1:4) .eq. '1aiu' .and. afil(ii)(1:5) .eq. '1r26a').or.  
     &( afil(1)(1:4) .eq. '1b3a' .and. afil(ii)(1:5) .eq. '1b3aa').or.  
     &( afil(1)(1:4) .eq. '1bkr' .and. afil(ii)(1:5) .eq. '1bkra').or.  
     &( afil(1)(1:4) .eq. '1c8c' .and. afil(ii)(1:5) .eq. '1xyia').or.  
     &( afil(1)(1:4) .eq. '1cei' .and. afil(ii)(1:5) .eq. '1fr2a').or.  
     &( afil(1)(1:4) .eq. '1cg5' .and. afil(ii)(1:5) .eq. '1irdb').or.  
!     & (afil(1)(1:4) .eq. '4ubp.pdb' .and. afil(ii) .eq. '1ay7b').or.
     &( afil(1)(1:4) .eq. '1dhn' .and. afil(ii)(1:5) .eq. '1nbua').or.  
     &( afil(1)(1:4) .eq. '1e6i' .and. afil(ii)(1:5) .eq. '1e6ia').or.  
     &( afil(1)(1:4) .eq. '1enh' .and. afil(ii)(1:5) .eq. '2hdda').or.  
     &( afil(1)(1:4) .eq. '1ew4' .and. afil(ii)(1:5) .eq. '1ew4a').or.  
     &( afil(1)(1:4) .eq. '1eyv' .and. afil(ii)(1:5) .eq. '1tzva').or.  
     &( afil(1)(1:4) .eq. '1fkb' .and. afil(ii)(1:5) .eq. '1jvwa').or.  
     &( afil(1)(1:4) .eq. '1fkb' .and. afil(ii)(1:5) .eq. '2awga').or.  
     &( afil(1)(1:4) .eq. '1fkb' .and. afil(ii)(1:5) .eq. '2f4ea').or.  
     &( afil(1)(1:4) .eq. '1hz6' .and. afil(ii)(1:5) .eq. '1hz6a').or.  
     &( afil(1)(1:4) .eq. '1lis' .and. afil(ii)(1:5) .eq. '2lisa').or.  
     &( afil(1)(1:4) .eq. '1lou' .and. afil(ii)(1:5) .eq. '1cqma').or.  
     &( afil(1)(1:4) .eq. '1nps' .and. afil(ii)(1:5) .eq. '1npsa').or.  
     &( afil(1)(1:4) .eq. '1ubi' .and. afil(ii)(1:5) .eq. '1yj1a').or.  
     &( afil(1)(1:4) .eq. '1ugh' .and. afil(ii)(1:5) .eq. '1ugia').or.  
     &( afil(1)(1:4) .eq. '1urn' .and. afil(ii)(1:5) .eq. '1nu4a').or.  
     &( afil(1)(1:4) .eq. '2acy' .and. afil(ii)(1:5) .eq. '1gxua').or.  
     &( afil(1)(1:4) .eq. '2chf' .and. afil(ii)(1:5) .eq. '1jbea').or.  
     &( afil(1)(1:4) .eq. '2ci2' .and. afil(ii)(1:5) .eq. '1vbwa').or.  
     &( afil(1)(1:4) .eq. '4ubp' .and. afil(ii)(1:5) .eq. '4ubpa').or.  
     &( afil(1)(1:4) .eq. '4pti' .and. afil(ii)(1:5) .eq. '1g6xa').or.  
     &( afil(1)(1:6).eq.'1hdd-C' .and. afil(ii)(1:5) .eq. '2hdda').or.  
     &( afil(1)(1:6).eq.'1bg8-A' .and. afil(ii)(1:5) .eq. '1dj8a').or.  
     &( afil(1)(1:4) .eq. '1bl0' .and. afil(ii)(1:5) .eq. '1bl0a').or.  
     &( afil(1)(1:4) .eq. '1jwe' .and. afil(ii)(1:5) .eq. '1b79a').or.  
     &( afil(1)(1:4) .eq. 'smd3' .and. afil(ii)(1:5) .eq. '1d3ba').or.  
     &( afil(1)(1:4) .eq. '1beo' .and. afil(ii)(1:5) .eq. '1lria').or.  
     &( afil(1)(1:4) .eq. '1nkl' .and. afil(ii)(1:5) .eq. '1l9la').or.  
     &( afil(1)(1:4) .eq. '1dtk' .and. afil(ii)(1:5) .eq. '1g6xa').or.  
     &( afil(1)(1:4) .eq. '2ovo' .and. afil(ii)(1:5) .eq. '1r0ri').or.  
     &( afil(1)(1:4) .eq. '1bbh' .and. afil(ii)(1:5) .eq. '1e85a').or.  
     &( afil(1)(1:4) .eq. '1bbh' .and. afil(ii)(1:5) .eq. '1s5aa').or.  
     &( afil(1)(1:4) .eq. '1c2r' .and. afil(ii)(1:5) .eq. '2bh4x').or.  
     &( afil(1)(1:4) .eq. '1cau' .and. afil(ii)(1:5) .eq. '1dgwx').or.  
     &( afil(1)(1:4) .eq. '1cau' .and. afil(ii)(1:5) .eq. '1dgwy').or.  
     &( afil(1)(1:4) .eq. '1cew' .and. afil(ii)(1:5) .eq. '1cewi').or.  
     &( afil(1)(1:4) .eq. '1dxt' .and. afil(ii)(1:5) .eq. '1irdb').or.  
     &( afil(1)(1:4) .eq. '1eaf' .and. afil(ii)(1:5) .eq. '1scza').or.  
     &( afil(1)(1:4) .eq. '1gky' .and. afil(ii)(1:5) .eq. '2an9a').or.  
     &( afil(1)(1:4) .eq. '1lga' .and. afil(ii)(1:5) .eq. '2boqa').or.  
     &( afil(1)(1:4) .eq. '2afn' .and. afil(ii)(1:5) .eq. '2bw4a').or.  
     &( afil(1)(1:4) .eq. '2fbj' .and. afil(ii)(1:5) .eq. '1mexl').or.  
     &( afil(1)(1:4) .eq. '2mta' .and. afil(ii)(1:5) .eq. '2c8sa').or.  
     &( afil(1)(1:4) .eq. '4sbv' .and. afil(ii)(1:5) .eq. '1c8na').or.  
     &( afil(1)(1:5).eq.'1b4bA' .and. afil(ii)(1:5) .eq. '1xxaa') .or.  
     &( afil(1)(1:5).eq. '1b72A' .and. afil(ii)(1:5) .eq. '2hdda').or.  
     &( afil(1)(1:5).eq. '1cewl' .and. afil(ii)(1:5) .eq. '1cewi').or.  
     &( afil(1)(1:5).eq. '1cqkA' .and. afil(ii)(1:5) .eq. '1oqoa').or.  
     &( afil(1)(1:5).eq. '1cy5A' .and. afil(ii)(1:5) .eq. '1cy5a').or.  
     &( afil(1)(1:5).eq. '1di2A' .and. afil(ii)(1:5) .eq. '1di2a').or.  
     &( afil(1)(1:5).eq. '1dtjA' .and. afil(ii)(1:5) .eq. '1zzka').or.  
     &( afil(1)(1:5).eq. '1g1cA' .and. afil(ii)(1:5) .eq. '1u2ha').or.  
     &( afil(1)(1:5).eq. '1gnuA' .and. afil(ii)(1:5) .eq. '1gnua').or.  
     &( afil(1)(1:5).eq. '1gpt_' .and. afil(ii)(1:5) .eq. '4sgbi').or.  
     &( afil(1)(1:5).eq. '1gyvA' .and. afil(ii)(1:5) .eq. '2a7ba').or.  
     &( afil(1)(1:5).eq. '1hbkA' .and. afil(ii)(1:5) .eq. '1hbka').or.  
     &( afil(1)(1:5).eq. '1jnuA' .and. afil(ii)(1:5) .eq. '1n91a').or.  
     &( afil(1)(1:5).eq. '1mn8A' .and. afil(ii)(1:5) .eq. '1mn8a').or.  
     &( afil(1)(1:5).eq. '1no5A' .and. afil(ii)(1:5) .eq. '1no5a').or.  
     &( afil(1)(1:5).eq. '1npsA' .and. afil(ii)(1:5) .eq. '1npsa').or.  
     &( afil(1)(1:6).eq. '1o2fB_' .and. afil(ii)(1:5).eq. '1svfa').or.  
     &( afil(1)(1:5).eq. '1orgA' .and. afil(ii)(1:5) .eq. '1ow4a')   
     &  )       
     &  ) then
             goto 3454     
          endif
!          print*, 'Normal'

          do i=1,20000
              xp(i)=0.0
              yp(i)=0.0
              zp(i)=0.0
          enddo
          do i=1,10000
              restyp(i) = 0
              resdep(i) = 0.0
          enddo
          do i=1,10000
              xp_n(i)=0.0
              yp_n(i)=0.0
              zp_n(i)=0.0
              xp_ca(i)=0.0
              yp_ca(i)=0.0
              zp_ca(i)=0.0
              xp_cb(i)=0.0
              yp_cb(i)=0.0
              zp_cb(i)=0.0
              xp_c(i)=0.0
              yp_c(i)=0.0
              zp_c(i)=0.0
              xp_o(i)=0.0
              yp_o(i)=0.0
              zp_o(i)=0.0
          enddo

            filename=afil(ii)
!            write(*,*) ii,' ',afil(ii)
         
          if(ii.eq.1) then

      open(11,file=fbase(1:ips)//'/2640total/'//file0(1:ips2),
     &  status='old')
          do j=1,10
             if(filename(j:j).eq.'.') idot=j-1
          enddo

      open(12,file=fbase(1:ips)//'/2640dep/'//file0(1:ips2)
     & //'_res.dep',status='old')
!      open(12,file=filename(1:5)//'_res.dep',status='old')
          else
      open(11,file=fbase(1:ips)//'/2640total/'//filename(1:5),
     & status='old')
      open(12,file=fbase(1:ips)//'/2640dep/'//filename(1:5)//'_res.dep',
     &             status='old')
          endif


!      print*, 'point 0'
       ime=-1 
       id=0
       id2=0
       id3=0
       id4=0
       id5=0
       id6=0
       numid=0
        numhl=0
        numst=0
        ires=0
        ares0='     '
        ares1='     '
        ares_n='     '
        ares_ca='     '
        ares_cb='     '
        ares_c='     '
        ares_o='     '
        ares_n_before='     '
        ares_n_before2='     '
        ares_n_before3='     '
        ares_ca_before='     '
        ares_ca_before2='     '
        ares_cb_before='     '
        ares_cb_before2='     '
        ares_c_before='     '
        ares_c_before2='     '
        ares_o_before='     '
        ares_o_before2='     '
        gly_flag = .FALSE.
             
15      read(11,'(a80)',err=10,end=10) whole1(id)
       if(whole1(id)(1:3).eq.'TER') goto 10

       if(whole1(id)(1:4).ne.'ATOM')   goto 15

       if(.not.(
     &(( whole1(id)(14:15) .eq. 'N ').or.
     &( whole1(id)(14:15) .eq. 'CA').or.
     &( whole1(id)(14:15) .eq. 'C ').or.
     &( whole1(id)(14:15) .eq. 'O ').or.
     &( whole1(id)(14:15) .eq. 'CB')) .and.
     &  ((whole1(id)(17:17) .eq. ' ') .or. (whole1(id)(17:17)
     &  .eq. 'A'))
     &         ) ) goto 15


      read(whole1(id)(18:21),'(a4)') resname 
!      READ (whole1(id)(31:54),'(3F8.3)') xptmp, yptmp,zptmp
      READ (whole1(id)(23:27),'(a5)') ares1

!       print*, 'resname',resname, xptmp,yptmp,zptmp, ares1
       if( whole1(id)(14:15) .eq. 'N ' .and.
     &  ((whole1(id)(17:17) .eq. ' ') .or. (whole1(id)(17:17) 
     &   .eq. 'A'))) then
!         print*, 'point 0'
        id2=id2+1
         ares_n_before2 = ares_n_before
         ares_n_before = ares_n
        READ (whole1(id)(31:54),'(3F8.3)') xp_n(id2),yp_n(id2),zp_n(id2)
        READ (whole1(id)(23:27),'(a5)') ares_n

        if( ares_n .ne. ares_n_before) then
          id=id+1
          xp(id) = xp_n(id2)
          yp(id) = yp_n(id2)
          zp(id) = zp_n(id2)

!          print*,'id',id,'N',xp(id),yp(id),zp(id)
        endif

       goto 5
       endif


       if( whole1(id)(14:15) .eq. 'CA' .and.
     &  ((whole1(id)(17:17) .eq. ' ') .or. (whole1(id)(17:17) 
     &   .eq. 'A'))) then
          id3=id3+1
!          ires=ires+1
          READ (whole1(id)(31:54),'(3F8.3)') xp_ca(id3),yp_ca(id3),
     &                                       zp_ca(id3)
          READ (whole1(id)(23:27),'(a5)') ares_ca
           do j1=1,20
            if(resn(j1).eq.resname(1:3)) then
              restyp(id3)=j1
            endif
           enddo

        if( (ares_ca .ne. ares_ca_before)
     &   .and. (ares_n .eq. ares_ca)
     &   ) then
          id=id+1
          xp(id) = xp_ca(id3)
          yp(id) = yp_ca(id3)
          zp(id) = zp_ca(id3)
!          print*,'id',id,'CA',xp(id),yp(id),zp(id)
          id=id+1
          xp(id) = xp_ca(id3)
          yp(id) = yp_ca(id3)
          zp(id) = xp_ca(id3)
!          id=id+1
!          xp(id) = xp_ca(id3)
!          yp(id) = yp_ca(id3)
!          zp(id) = xp_ca(id3)
!          id=id+1
!          xp(id) = xp_ca(id3)
!          yp(id) = yp_ca(id3)
!          zp(id) = xp_ca(id3)
        endif

       goto 5
       endif


       if( whole1(id)(14:15) .eq. 'C ' .and.
     &  ((whole1(id)(17:17) .eq. ' ') .or. (whole1(id)(17:17) 
     &  .eq. 'A'))) then
!         print*, 'point 0'
          id4=id4+1
          READ (whole1(id)(31:54),'(3F8.3)') xp_c(id4),yp_c(id4),
     &                                       zp_c(id4)
          READ (whole1(id)(23:27),'(a5)') ares_c

        if(ares_n .eq. ares_c .and. ares_c .ne. ares_c_before) then
          id=id+1
          xp(id) = xp_c(id4)
          yp(id) = yp_c(id4)
          zp(id) = zp_c(id4)
!          print*,'id',id,'C',xp(id),yp(id),zp(id)
        endif

       goto 5
       endif


       if( whole1(id)(14:15) .eq. 'O ' .and.
     &  ((whole1(id)(17:17) .eq. ' ') .or. (whole1(id)(17:17) 
     &  .eq. 'A'))) then
!         print*, 'point 0'
          id5=id5+1
          READ (whole1(id)(31:54),'(3F8.3)') xp_o(id5),yp_o(id5),
     &                                       zp_o(id5)
          READ (whole1(id)(23:27),'(a5)') ares_o

        if(ares_n .eq. ares_o .and. ares_o .ne. ares_o_before) then
          id=id+1
          xp(id) = xp_o(id5)
          yp(id) = yp_o(id5)
          zp(id) = zp_o(id5)
!          print*,'id',id,'O',xp(id),yp(id),zp(id)
        endif

       goto 5
       endif


       if( whole1(id)(14:15) .eq. 'CB' .and.
     &  ((whole1(id)(17:17) .eq. ' ') .or. (whole1(id)(17:17) 
     &  .eq. 'A'))) then
          id6=id6+1
          READ (whole1(id)(31:54),'(3F8.3)') xp_cb(id6),yp_cb(id6),
     &                                       zp_cb(id6)
          READ (whole1(id)(23:27),'(a5)') ares_cb

        if( ares_n .eq. ares_cb 
     &       .and. ares_cb .ne. ares_cb_before) then
          id=id-2
          xp(id) = xp_cb(id6)
          yp(id) = yp_cb(id6)
          zp(id) = zp_cb(id6)
!          print*,'id',id,'CB',xp(id),yp(id),zp(id)
          id=id+2
        endif

       goto 5
       endif


 5       if(
     &       ( (ares_n .ne. ares_ca) .and. (ares_n .ne. ares_c) .and.
     &         (ares_n .ne. ares_o) .and.
!     &         (ares_ca .eq. ares_c) .and. (ares_ca .eq. ares_o)
     &         (id .ge. 5) )
     &     ) then

        if( (ares_ca .ne. ares_cb) .and.
     &      (ares_ca .eq. ares_c) .and.
     &      (ares_ca .eq. ares_o)
     &    ) then

!      print*, 'point 1'
          cb(1)=0.0
          cb(2)=0.0
          cb(3)=0.0

          x2(1)=xp_ca(id3)
          x2(2)=yp_ca(id3)
          x2(3)=zp_ca(id3)

          x3(1)=xp_n(id2-1)
          x3(2)=yp_n(id2-1)
          x3(3)=zp_n(id2-1)
          x4(1)=xp_c(id4)
          x4(2)=yp_c(id4)
          x4(3)=zp_c(id4)

!         print*,'ii',ii, 'id', id, x2(1), x2(2), x2(3)
        call  xyzatm (cb,x2,bond,x4,angle1,x3,angle2,0)
!         print*,cb(1), cb(2), cb(3)

           xp(id-3)=cb(1)
           yp(id-3)=cb(2)
           zp(id-3)=cb(3)

          id6=id6+1
          xp_cb(id6) = cb(1)
          yp_cb(id6) = cb(2)
          zp_cb(id6) = cb(3)
       endif


       ares_n_before3 = ares_n_before2
       ares_ca_before2 = ares_ca_before
       ares_cb_before2 = ares_cb_before
       ares_c_before2 = ares_c_before
       ares_o_before2 = ares_o_before

       ares_n_before2 = ares_n_before
       ares_ca_before = ares_ca
       ares_cb_before = ares_cb
       ares_c_before = ares_c
       ares_o_before = ares_o
        endif



        if(ares1.ne.ares0) then
           ires=ires+1
           ares0=ares1
!           restyp(ires)=-1
        endif

       natom=id3
       nres=id3
       rnum=ires

      goto 15

 10   continue

       if( (ares_ca .ne. ares_cb)  .and.
     &      (ares_ca .eq. ares_c) .and.
     &      (ares_ca .eq. ares_o)
     &    ) then

          cb(1)=0.0
          cb(2)=0.0
          cb(3)=0.0

          x2(1)=xp_ca(id3)
          x2(2)=yp_ca(id3)
          x2(3)=zp_ca(id3)

          x3(1)=xp_n(id2)
          x3(2)=yp_n(id2)
          x3(3)=zp_n(id2)

          x4(1)=xp_c(id4)
          x4(2)=yp_c(id4)
          x4(3)=zp_c(id4)


        call  xyzatm (cb,x2,bond,x4,angle1,x3,angle2,0)

           xp(id-2)=cb(1)
           yp(id-2)=cb(2)
           zp(id-2)=cb(3)
        endif

       close(11)

           do i=1,rnum
              read(12,*,err=20,end=20) j,tmp
              resdep(i)=exp(-tmp*dd)
           enddo
 20        continue
           close(12)


!      print*,'id2',id2,'id3',id3,'id4',id4,'id5',id5,'id6',id6,'id',id
!       print*, 'rnum', rnum,'natom',natom
        if(ii.eq.1) then
         rnumme=rnum
         do i=1,natom*5
         xp0(i)=xp(i)
         yp0(i)=yp(i)
         zp0(i)=zp(i)
         enddo
         do i=1,natom
         mydep(i)=resdep(i)
         restyp1(i)=restyp(i)
         enddo
         write(*,'(70a)') (resmap(restyp1(i)),i=1,rnumme)

        else 


        do isf=1,rnumme-frgwind+1
           do j=1,frgwind
           lr0=j+isf-1
           nn=15*(j-1)+1
           xx(nn)=xp0(5*(lr0-1)+1)
           xx(nn+1)=yp0(5*(lr0-1)+1)
           xx(nn+2)=zp0(5*(lr0-1)+1)
           xx(nn+3)=xp0(5*(lr0-1)+2)
           xx(nn+4)=yp0(5*(lr0-1)+2)
           xx(nn+5)=zp0(5*(lr0-1)+2)
           xx(nn+6)=xp0(5*(lr0-1)+3)
           xx(nn+7)=yp0(5*(lr0-1)+3)
           xx(nn+8)=zp0(5*(lr0-1)+3)
           xx(nn+9)=xp0(5*(lr0-1)+4)
           xx(nn+10)=yp0(5*(lr0-1)+4)
           xx(nn+11)=zp0(5*(lr0-1)+4)
           xx(nn+12)=xp0(5*(lr0-1)+5)
           xx(nn+13)=yp0(5*(lr0-1)+5)
           xx(nn+14)=zp0(5*(lr0-1)+5)
           enddo
        do i=1,rnum-frgwind+1
         disp=0
         dispa=0.0
        do j=1,frgwind 
           lr=j+i-1
           nn=15*(j-1)+1
           yy(nn)=xp(5*(lr-1)+1)
           yy(nn+1)=yp(5*(lr-1)+1)
           yy(nn+2)=zp(5*(lr-1)+1)
           yy(nn+3)=xp(5*(lr-1)+2)
           yy(nn+4)=yp(5*(lr-1)+2)
           yy(nn+5)=zp(5*(lr-1)+2)
           yy(nn+6)=xp(5*(lr-1)+3)
           yy(nn+7)=yp(5*(lr-1)+3)
           yy(nn+8)=zp(5*(lr-1)+3)
           yy(nn+9)=xp(5*(lr-1)+4)
           yy(nn+10)=yp(5*(lr-1)+4)
           yy(nn+11)=zp(5*(lr-1)+4)
           yy(nn+12)=xp(5*(lr-1)+5)
           yy(nn+13)=yp(5*(lr-1)+5)
           yy(nn+14)=zp(5*(lr-1)+5)
        enddo

!      print*, 'point 2'
        call fitsq(disp,xx,yy,frgwind*5,tt,bb,convg)
!      print*, 'point 2.1'
           dispa=0.0
           do j=1,frgwind
             lr0=j+isf-1
             lr=j+i-1
             dispa=dispa+(resdep(lr) - mydep(lr0))**2
           enddo
          disp=disp+dispa*cov
          
!      print*, 'point 2.5'
        call pushe(disp,disp0,
     &    inb,isf,i,rnumme,frgwind,restyp, typhold, poshld,
     &    filename,hlc)

!      print*, 'point 3'

        enddo 
        enddo

        if(mod(ii,10) .eq. 0) then
       write(*,*) afil(ii),rnum,ii
        endif
         
        endif
       
 3454  continue

       enddo  ! ii


!         print*,'point 4'
          misf=rnumme-frgwind+1 

          write(24,*) misf,inb,frgwind,rnumme
          do isf1=1,misf

          do i=1,inb
           disp1(i,isf1)=sqrt(disp0(i,isf1))
           write(24,*) isf1,i,disp1(i,isf1),
     &    (typhold(i,k,isf1),k=1,frgwind),
     &    '  ',hlc(i,isf1),' ',poshld(i,isf1)
c           write(24,*) mydep(isf1),(hold(i,k,isf1),k=1,frgwind)
          enddo            
          enddo

!         print*,'point 5'
          close(24)

! 500    continue
c     now construct profile and get the most probable sequence
!        write(91,*) rnumme,nfil
c 3455   continue
        idt2_60=0.0
        idt2_65=0.0
        idt2_80=0.0

        print*, 'now calling secondpro'
        radius = 6.0
        call secondpro(fbase, file0, ips, ips2,nfil, afil, fa2_60,
     &       idt2_60, radius)

        print*, 'now calling thirdpro'
        radius = 8.0
        call secondpro(fbase, file0, ips, ips2,nfil, afil,fa2_80,
     &       idt2_80, radius)

        print*, 'idt2 60, 80', idt2_60, idt2_80       
        if(idt2_60 > idt2_80) then
             idt2 = idt2_60 
           do i=1,rnumme
             do j=1,20
               fa2(j,i)=fa2_60(j,i)
             enddo
           enddo
        else if(idt2_80 >= idt2_60 ) then
             idt2 = idt2_80
           do i=1,rnumme
             do j=1,20
               fa2(j,i)=fa2_80(j,i)
             enddo
           enddo
        endif 

        print*, 'idt2 ', idt2       


        do i=1,rnumme

            do j=1,20
               fa(j,i)=0.0
               fa3(j,i)=0.0
            enddo

         do k=1,frgwind
             l=i-k+1
             if(l.gt.0.and.l.le.misf) then
                  do n=1,inb
                  j=typhold(n,k,l)
!                 temp=real(floor(real(n)/5.0)+1)
!                 fa(j,i)=fa(j,i)+1.0*(1.0/temp)

                  fa(j,i)=fa(j,i)+1.0
                  enddo
             endif
         enddo

             tmp=0.0
             tmp1=-1.0
             imax1(i)=-1
            do j=1,20
              if(fa(j,i).ge.tmp1) then
                 tmp1=fa(j,i)
                imax1(i)=j
              endif
              tmp=tmp+fa(j,i)
            enddo
            tmp=1./tmp
            do j=1,20
              fa(j,i)=fa(j,i)*tmp
            enddo


             tmp=0.0
             tmp1=-1.0
             imax2(i)=-1
            do j=1,20
              if(fa2(j,i).ge.tmp1) then
                 tmp1=fa2(j,i)
                 imax2(i)=j
              endif
!              tmp=tmp+fa2(j,i)
            enddo

!             write(*,*) 'fa2',(fa2(j,i),j=1,20)
!            tmp=1./tmp
!            do j=1,20
!              fa2(j,i)=fa2(j,i)*tmp
!            enddo


            
            do j=1,20
!!        fa3(j,i)=alp*fa(j,i)+(1-alp)*fa2(j,i) + 1.25*fa(j,i)*fa2(j,i)
         fa3(j,i)=alp*fa(j,i)+(1.0-alp)*fa2(j,i)
            enddo

!             write(*,*) 'fa3',(fa3(j,i),j=1,20)


             tmp=0.0
             tmp1=-10000.0
             imax3(i)=-1
            do j=1,20
              if(fa3(j,i).ge.tmp1) then
                 tmp1=fa3(j,i)
                 imax3(i)=j
              endif
              tmp=tmp+fa3(j,i)
            enddo
            tmp=1./tmp


            it=restyp1(i)
         write(91,100) i,' ',resmap(it),' ',it,
     &   (fa3(j,i)*tmp,j=1,20)
 100   format (I6, A, A, A, I3, 2X, 20F10.5)

!         print*,'point 4'
        enddo    

          close(91)

!         print*,'point 5'
c     write out seq:
         idt1 = 0.0
         idt2 = 0.0
         idt3 = 0.0
         idt1_gcp = 0.0
         idt2_gcp = 0.0
         idt3_gcp = 0.0
         idt1_nongcp = 0.0
         idt2_nongcp = 0.0
         idt3_nongcp = 0.0
         write(*,'(a1,a5,i5)') '>',file0,rnumme
         write(*,'(70a)') (resmap(restyp1(i)),i=1,rnumme)

         print*,'imax1'
         write(*,'(70a)') (resmap(imax1(i)),i=1,rnumme)
         print*,'imax2_after'
         write(*,'(70a)') (resmap(imax2(i)),i=1,rnumme)
         print*,'imax3'
         write(*,'(70a)') (resmap(imax3(i)),i=1,rnumme)


         rnumme_gcp=0.0
         rnumme_nongcp=0.0
         do i=1, rnumme
           if( (resmap(restyp1(i)) .eq. 'G') .or.
     &          (resmap(restyp1(i)) .eq. 'C') .or.
     &          (resmap(restyp1(i)) .eq. 'P') ) then
              rnumme_gcp=rnumme_gcp + 1.0
              if(restyp1(i) .eq. imax1(i) ) then
                idt1=idt1+1.0
                idt1_gcp=idt1_gcp + 1.0
              endif
           else
              rnumme_nongcp=rnumme_nongcp + 1.0
              if(restyp1(i) .eq. imax1(i) ) then
                idt1=idt1+1.0
                idt1_nongcp=idt1_nongcp + 1.0
              endif
           endif
         enddo

         write(*,*) '      '
         idt1 = idt1*100.0/((rnumme)*1.0)
         idt1_gcp = idt1_gcp*100.0/rnumme_gcp
         idt1_nongcp = idt1_nongcp*100.0/rnumme_nongcp
         write(*,777) file0,'  idt1 is ',idt1
         write(*,*) '      '



         rnumme_gcp=0.0
         rnumme_nongcp=0.0
         do i=1, rnumme
           if( (resmap(restyp1(i)) .eq. 'G') .or.
     &          (resmap(restyp1(i)) .eq. 'C') .or.
     &          (resmap(restyp1(i)) .eq. 'P') ) then
              rnumme_gcp=rnumme_gcp + 1.0
              if(restyp1(i) .eq. imax2(i) ) then
                idt2=idt2+1.0
                idt2_gcp=idt2_gcp + 1.0
              endif
           else
              rnumme_nongcp=rnumme_nongcp + 1.0
              if(restyp1(i) .eq. imax2(i) ) then
                idt2=idt2+1.0
                idt2_nongcp=idt2_nongcp + 1.0
              endif
           endif
         enddo

         write(*,*) '      '
         idt2 = idt2*100.0/((rnumme)*1.0)
         idt2_gcp = idt2_gcp*100.0/rnumme_gcp
         idt2_nongcp = idt2_nongcp*100.0/rnumme_nongcp
         write(*,777) file0,'  idt2 is ', idt2
         write(*,*) '      '




         rnumme_gcp=0.0
         rnumme_nongcp=0.0
         do i=1, rnumme
           if( (resmap(restyp1(i)) .eq. 'G') .or.
     &          (resmap(restyp1(i)) .eq. 'C') .or.
     &          (resmap(restyp1(i)) .eq. 'P') ) then
              rnumme_gcp=rnumme_gcp + 1.0
              if(restyp1(i) .eq. imax3(i) ) then
                idt3=idt3+1.0
                idt3_gcp=idt3_gcp + 1.0
              endif
           else
              rnumme_nongcp=rnumme_nongcp + 1.0
              if(restyp1(i) .eq. imax3(i) ) then
                idt3=idt3+1.0
                idt3_nongcp=idt3_nongcp + 1.0
              endif
           endif
         enddo


         write(*,*) '      '
         idt3 = idt3*100.0/((rnumme)*1.0)
         idt3_gcp = idt3_gcp*100.0/rnumme_gcp
         idt3_nongcp = idt3_nongcp*100.0/rnumme_nongcp
         write(*,777) file0,'  idt3 is ', idt3
         write(*,*) '      '


         write(*,*) file0,'  idt1+idt2 is ',(idt1+idt2) 
         write(*,*) '      '



 777   format (A20, A10, F6.3)

       stop
      end


      subroutine msort(fa,orderml,ii,numt)
      integer fa(35,1070),orderml(50),ii,numt(1070)

      logical corre2(50),yes
      real*8 mintmp
      integer minhld1,k,k1,k2

          do k1=1,numt(ii)
              orderml(k1)=k1
          enddo

          do k=1,numt(ii)
             corre2(k)= .FALSE.
          enddo
          temp = 0
          do k1=1,numt(ii)
                mintmp = 10000
                minhld1 = 0
                yes = .FALSE.
           do k2 = 1,numt(ii)
              if ( (fa(k1,ii) .gt. 0 ) .and. (fa(k2,ii) .gt.
     &            0) .and. ( fa(k2,ii) .lt. mintmp ) .and.
     &            ( corre2(k2) .eqv. .FALSE.) )  then
                 mintmp = fa(k2,ii)
                 minhld1 = k2
                 yes = .TRUE.
              endif
           enddo
                 if ( yes .eqv. .TRUE.) then
                    orderml(k1) = minhld1
                    corre2(minhld1) = .TRUE.
                 endif
          enddo

         return
         end


      subroutine msort2(fa,orderml,numl)
      integer fa(35),orderml(50),numl

      logical corre2(50),yes
      real*8 mintmp
      integer minhld1,k,k1,k2

          do k1=1,numl
              orderml(k1)=k1
          enddo

          do k=1,numl
             corre2(k)= .FALSE.
          enddo
          temp = 0
          do k1=1,numl
                mintmp = 10000
                minhld1 = 0
                yes = .FALSE.
           do k2 = 1,numl
              if ( (fa(k1) .gt. 0 ) .and. (fa(k2) .gt.
     &            0) .and. ( fa(k2) .lt. mintmp ) .and.
     &            ( corre2(k2) .eqv. .FALSE.) )  then
                 mintmp = fa(k2)
                 minhld1 = k2
                 yes = .TRUE.
              endif
           enddo
                 if ( yes .eqv. .TRUE.) then
                    orderml(k1) = minhld1
                    corre2(minhld1) = .TRUE.
                 endif
          enddo

         return
         end

      
      subroutine mysort(fa,orderml,ii)
      real*8 fa(20,10000)
      integer orderml(20),ii

      logical corre2(20),yes
      real*8 maxtmp
      integer maxhld1,k,k1,k2

          do k1=1,20
              orderml(k1)=k1
          enddo

          do k=1,20
             corre2(k)= .FALSE.
          enddo
          temp = 0
          do k1=1,20
                maxtmp = -1.0
                maxhld1 = 0
                yes = .FALSE.
           do k2 = 1,20
              if ( (fa(k1,ii) .gt. 0.0 ) .and. (fa(k2,ii) .gt.
     &            0.0) .and. ( fa(k2,ii) .gt. maxtmp ) .and.
     &            ( corre2(k2) .eqv. .FALSE.) )  then
                 maxtmp = fa(k2,ii)
                 maxhld1 = k2
                 yes = .TRUE.
              endif
           enddo
                 if ( yes .eqv. .TRUE.) then
                    orderml(k1) = maxhld1
                    corre2(maxhld1) = .TRUE.
                 endif
          enddo

         return
         end



          subroutine pushe(disp,disp0,
     &    inb, isf,i1,rnumme,frgwind,restyp, typhold0, poshld,
     &    filename,hlc)
          real*8 disp, disp0(25,10000)
          integer inb, isf, i1, rnumme,frgwind,restyp(10000),
     &    typhold0(25,9,10000),poshld(25,10000)
          character filename*40,hlc(25,10000)*10

          
          real emax 
          integer me, i, j, jmax

           emax=-1.d8
           me=isf
        do i=1,inb
           if(disp0(i,me).gt.emax) then
             jmax=i
             emax=disp0(i,me)
           endif
        enddo


          if(disp.lt.emax) then 
             do j=1,frgwind
                i=i1+j-1
             typhold0(jmax,j,me)=restyp(i)
             enddo
             disp0(jmax,me)=disp
             hlc(jmax,me)=filename
             poshld(jmax,me)=i1
             
!             print*, disp
!             print*, (typhold0(jmax,j,me),j=1,frgwind)
          endif

          return
          end




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine fitsq2(rms,x,y,tc,lc,nn,t,b,non_conv)
      implicit real*8 (a-h,o-z)
c      include 'COMMON.IOUNITS'
c  x and y are the vectors of coordinates (dimensioned (3,n)) of the two
c  structures to be superimposed.  nn is 3*n, where n is the number of  
c  points.   t and b are respectively the translation vector and the    
c  rotation matrix that transforms the second set of coordinates to the 
c  frame of the first set.                                              
c  eta =  machine-specific variable                                     
                                                                        
      real*8 rms, x(30),y(30),tc(3), lc(3)               
      integer nn
      real*8 t(3),b(3,3)
      logical non_conv

      dimension q(3,3),r(3,3),v(3),xav(3),yav(3),e(3),c(3,3)     
      eta = z00100000                                                   
c     small=25.0*rmdcon(3)                                              
c     small=25.0*eta                                                    
c     small=25.0*10.e-10                                                
c the following is a very lenient value for 'small'                     
      small = 0.0001D0                                                  
      non_conv=.false.
      fn=nn                                                             
      do 10 i=1,3                                                       
      xav(i)=0.0D0                                                      
      yav(i)=0.0D0                                                      
      do 10 j=1,3                                                       
   10 b(j,i)=0.0D0                                                      
      nc=0                                                              
c                                                                       
!      do 30 n=1,nn                                                      
!      do 20 i=1,3                                                       
c      write(iout,*)'x = ',x(nc+i),'  y = ',y(nc+i)                           
!      xav(i)=xav(i)+x(nc+i)/fn                                          
!   20 yav(i)=yav(i)+y(nc+i)/fn                                          
!   30 nc=nc+3                                                           
                                                                       
      do i=1,3                                                       
      xav(i)=tc(i)                                          
      yav(i)=lc(i)                                         
      enddo

      do i=1,3
        t(i)=yav(i)-xav(i)
      enddo
!        write(*,*) t(1), t(2), t(3)
!        write(*,*) '    ' 

      rms=0.0d0
      do n=1,nn
        do i=1,3
          rms=rms+(y(3*(n-1)+i)-x(3*(n-1)+i)-t(i))**2
        enddo
      enddo
      rms=dabs(rms/fn)

c     write(iout,*)'xav = ',(xav(j),j=1,3)                                    
c     write(iout,*)'yav = ',(yav(j),j=1,3)                                    
c     write(iout,*)'t   = ',(t(j),j=1,3)
c     write(iout,*)'rms=',rms
      if (rms.lt.small) return
                                                                        
                                                                        
      nc=0                                                              
      rms=0.0D0                                                         
      do 50 n=1,nn                                                      
      do 40 i=1,3                                                       
      rms=rms+((x(nc+i)-xav(i))**2+(y(nc+i)-yav(i))**2)/fn              
      do 40 j=1,3                                                       
      b(j,i)=b(j,i)+(x(nc+i)-xav(i))*(y(nc+j)-yav(j))/fn                
   40 c(j,i)=b(j,i)                                                     
   50 nc=nc+3                                                           
      call sivade(b,q,r,d,non_conv)
      !  write(*,*) t(1), t(2), t(3)
      !  write(*,*) '    ' 
      sn3=dsign(1.0d0,d)                                                   
      do 120 i=1,3                                                      
      do 120 j=1,3                                                      
  120 b(j,i)=-q(j,1)*r(i,1)-q(j,2)*r(i,2)-sn3*q(j,3)*r(i,3)             
      call mvvad(b,xav,yav,t)                                           
      !  write(*,*) t(1), t(2), t(3)
      !  write(*,*) '    ' 
      do 130 i=1,3                                                      
      do 130 j=1,3                                                      
      rms=rms+2.0*c(j,i)*b(j,i)                                         
  130 b(j,i)=-b(j,i)                                                    
      if (dabs(rms).gt.small) go to 140                                  
*     write (6,301)                                                     
      return                                                            
  140 if (rms.gt.0.0d0) go to 150                                         
c     write (iout,303) rms                                                 
      rms=0.0d0
*     stop                                                              
c 150 write (iout,302) dsqrt(rms)                                           
  150 continue
      !  write(*,*) t(1), t(2), t(3)
      !  write(*,*) '    ' 
      return                                                            
  301 format (5x,'rms deviation negligible')                            
  302 format (5x,'rms deviation ',f14.6)                                
  303 format (//,5x,'negative ms deviation - ',f14.6)                   
      end                                                               


      subroutine secondpro(fbase,file0, ips2, ips3, nfil,afil,fa2,
     & idt2_60,radius)
      implicit real*8(a-h,o-z)
      character fbase*80, file0*80
      integer ips2, ips3
      integer nfil
      character afil(10000)*80
      real*8 fa2(20,2000)
      real*8 idt2_60
      real*8 radius

      integer maxa, maxit, inb0, inb, inbb,ntyp, natyp, maxres
      real*8  cbox
      parameter(maxa=1270,cbox=18.856,maxit=100,inb0=25,inb=480,
     &          ntyp=20,natyp=8,maxres=1990)
      logical convg, corres1(inb, 37,1990),corres2(inb,37,1990),
     &        map(37)
      integer natom,natome,numdel,
     &        frgwind,m, isf,lhold(480,maxres),
     &        numlhold(150,1990),orderml(150,1990),
     &        restyp(1270),rnum,
     &        rnumme,restyp1(maxres),
     &        typhold(inb,37,maxres),imax(maxres),
     &        poshld(inb,maxres),numt(maxres),numl,
     &        distl(inb,37,maxres),jmax,
     &        i_count,nres,
     & tmp_resnumt(9),tmp_resnuml(9),resnumt(1990),resnuml(1270),
     & resnum_thold(37,1990),
     & resnum_lhold(480,37,1990),sel_resnuml(37),
     & resconta_thold(300,1990),l,l2,l3,n,i,j,k,id,id2,numid,o
     & ,idtemp,lr0,ii

      real*8 pi,pi3,rd,xxp(3750),xp(maxa),yyp(3750),yp(maxa),
     &       zzp(3750),zp(maxa),sphxt(37,maxres),
     &       sphyt(37,maxres), sphzt(37,maxres), 
     &       xpl(37),ypl(37),zpl(37),tmpt(3), tmpl(3), 
     &       tt(3),bb(3,3),fa(20,maxres),
     &   xx(30),txx(30),yy(30),tyy(30),xxp0(3750),xp0(maxa),
     &   yyp0(3750),yp0(maxa),zzp0(3750),zp0(maxa),
     &sphholdx(inb,37,maxres),sphholdy(inb,37,maxres),
     &sphholdz(inb,37,maxres),bbhold(inb,9,maxres), 
     &tthold(inb, 3, maxres),disp,disp0(480,maxres)
     &  ,disp1(480,maxres)
     &  ,centert(3,1990),centerl(3,1270),ctrlhold(480,3,1990)
     &  ,distml(150,1990), dist(37), disl(37), tmp0(37,3),tmp2(37,3) 
     &  ,xpmean,ypmean,zpmean,xptmp,yptmp,zptmp,temp,idt
     &  ,disttot(1990),mindist,resdep_lhold(480,37,1990)
     &  ,resdep(maxres),mydep(maxres)
     &  ,dispa


      character resname*4,resn(ntyp)*3,
     &          ares1*5,ares0*5,aresbefore*5,
c     &          hlc(480,maxres)*10,
     &          resmap(20)*1
      character whole*80, filename*40,whole1(maxa)*80,
     &      fil1*40,
     &       fil2*40, tmpl_name*40, dbf_name*40,
     &       pdbt*40, pdbl*40,resnl(37)*3, resnt(37)*3,
     &       xxresnt(9)*3, yyresnl(9)*3, cent_thold(maxres)*3,
     &       cent_tmp_l*3, cent_lhold(inb,maxres)*3

      data pi/3.1415926/

      data resn/'CYS','MET','PHE','ILE','LEU','VAL','TRP','TYR',
     &          'ALA','GLY','THR','SER','GLN','ASN','GLU','ASP',
     &          'HIS','ARG','LYS','PRO'/
   

      data resmap/'C','M','F','I','L','V','W','Y','A','G',
     &            'T','S','Q','N','E','D','H','R','K','P'/

c     
!        read(*, '(a80)')  fbase
        ips=index(fbase,' ')-1

!
!        i=1 
!
! 2001   read(*,*,end=2000) afil(i)
!        i=i+1
!        goto 2001
!
! 2000   continue
!        nfil=i-1

!         call getarg(1,fbase)
!         call getarg(2,file0)

!           read(*, '(a80)')  fbase
!         ips=index(fbase,' ')-1
!         fbase=trim(fbase)
!         file0=trim(file0)

!        open(11,file=fbase(1:ips)//'/lib-bin/testpdb.lst',status='old')

!        i=1
!        afil(1)=file0
!        i=i+1
c        open(unit=20,file='pdbfile.list',status='old')
! 2001   read(11,*,end=2000) afil(i)
!        i=i+1
!        goto 2001

! 2000   continue
!        nfil=i-1
c        close(20)


       convg = .FALSE.
       i_count = 0
       frgwind=9
       dd=1./2.8
       cov=10.0

          do i=1,inb
               do j=1,maxres
               disp0(i,j)=1.d8
               enddo
          enddo

          do i=1,480
            do k=1, maxres
               do j=1,37
                    typhold(i,j,k)=0
                    sphholdx(i,j,k)=0.0
                    sphholdy(i,j,k)=0.0
                    sphholdz(i,j,k)=0.0
                    distl(i,j,k) = 0
                    corres1(i,j,k) = .FALSE.
                    corres2(i,j,k) = .FALSE.
                    resnum_lhold(i,j,k) = 0
                    resdep_lhold(i,j,k) = 0.0
               enddo
               do j=1,9
                    bbhold(i,j,k)=0.0
               enddo
               do j=1,3
                    tthold(i,j,k)=0.0
                    ctrlhold(i,j,k)=0.0
               enddo
                    lhold(i,k)=0
                    cent_lhold(i,k)='   '
             enddo
          enddo

          do i=1,37
            do k=1,maxres
               sphxt(i,k)=0.0
               sphyt(i,k)=0.0
               sphzt(i,k)=0.0
               numlhold(i,k)=0
               orderml(i,k)=0
               distml(i,k)=0
            enddo
            map(i) = .FALSE.
            sel_resnuml(i) = 0
          enddo

          do i=1,300
            do k=1,maxres
               resconta_thold(i,k)=0
            enddo
          enddo

          do i=1,37
            do k=1,1990
               resnum_thold(i,k)=0
            enddo
          enddo

          do i=1,maxres
              numt(i)=0
              cent_thold(i) = '   '
              resnumt(i) = 0
             do k=1, 3
              centert(k,i)=0.0
             enddo
             resdep(i)=0.0
             mydep(i)=0.0
          enddo 

          do i=1,3750
              xxp0(i)=0.0
              yyp0(i)=0.0
              zzp0(i)=0.0
          enddo

      do ii=1,nfil

c         if(afil(1) .eq. '4ubp.pdb' .and. afil(ii) .eq. '1ay7b')then
c             print*, afil(ii),'Yes'
c         else 
c             print*, afil(ii),'No'
c         endif

 
          do i=1,3750
              xxp(i)=0.0
              yyp(i)=0.0
              zzp(i)=0.0
          enddo

          do i=1,1270
             do k=1, 3
              centerl(k,i)=0.0
             enddo
              resnuml(i)= 0
              restyp(i)= 0
              whole1(i)=' '
              xp(i)=0.0
              yp(i)=0.0
              zp(i)=0.0
          enddo 

         if ( (ii .ne. 1) .and. 
     &      (
     &( afil(1)(1:4) .eq. '1a19' .and. afil(ii)(1:5) .eq. '1ay7b').or.  
     &( afil(1)(1:4) .eq. '1aiu' .and. afil(ii)(1:5) .eq. '1gh2a').or.  
     &( afil(1)(1:4) .eq. '1aiu' .and. afil(ii)(1:5) .eq. '1r26a').or.  
     &( afil(1)(1:4) .eq. '1b3a' .and. afil(ii)(1:5) .eq. '1b3aa').or.  
     &( afil(1)(1:4) .eq. '1bkr' .and. afil(ii)(1:5) .eq. '1bkra').or.  
     &( afil(1)(1:4) .eq. '1c8c' .and. afil(ii)(1:5) .eq. '1xyia').or.  
     &( afil(1)(1:4) .eq. '1cei' .and. afil(ii)(1:5) .eq. '1fr2a').or.  
     &( afil(1)(1:4) .eq. '1cg5' .and. afil(ii)(1:5) .eq. '1irdb').or.  
!     & (afil(1)(1:4) .eq. '4ubp.pdb' .and. afil(ii) .eq. '1ay7b').or.
     &( afil(1)(1:4) .eq. '1dhn' .and. afil(ii)(1:5) .eq. '1nbua').or.  
     &( afil(1)(1:4) .eq. '1e6i' .and. afil(ii)(1:5) .eq. '1e6ia').or.  
     &( afil(1)(1:4) .eq. '1enh' .and. afil(ii)(1:5) .eq. '2hdda').or.  
     &( afil(1)(1:4) .eq. '1ew4' .and. afil(ii)(1:5) .eq. '1ew4a').or.  
     &( afil(1)(1:4) .eq. '1eyv' .and. afil(ii)(1:5) .eq. '1tzva').or.  
     &( afil(1)(1:4) .eq. '1fkb' .and. afil(ii)(1:5) .eq. '1jvwa').or.  
     &( afil(1)(1:4) .eq. '1fkb' .and. afil(ii)(1:5) .eq. '2awga').or.  
     &( afil(1)(1:4) .eq. '1fkb' .and. afil(ii)(1:5) .eq. '2f4ea').or.  
     &( afil(1)(1:4) .eq. '1hz6' .and. afil(ii)(1:5) .eq. '1hz6a').or.  
     &( afil(1)(1:4) .eq. '1lis' .and. afil(ii)(1:5) .eq. '2lisa').or.  
     &( afil(1)(1:4) .eq. '1lou' .and. afil(ii)(1:5) .eq. '1cqma').or.  
     &( afil(1)(1:4) .eq. '1nps' .and. afil(ii)(1:5) .eq. '1npsa').or.  
     &( afil(1)(1:4) .eq. '1ubi' .and. afil(ii)(1:5) .eq. '1yj1a').or.  
     &( afil(1)(1:4) .eq. '1ugh' .and. afil(ii)(1:5) .eq. '1ugia').or.  
     &( afil(1)(1:4) .eq. '1urn' .and. afil(ii)(1:5) .eq. '1nu4a').or.  
     &( afil(1)(1:4) .eq. '2acy' .and. afil(ii)(1:5) .eq. '1gxua').or.  
     &( afil(1)(1:4) .eq. '2chf' .and. afil(ii)(1:5) .eq. '1jbea').or.  
     &( afil(1)(1:4) .eq. '2ci2' .and. afil(ii)(1:5) .eq. '1vbwa').or.  
     &( afil(1)(1:4) .eq. '4ubp' .and. afil(ii)(1:5) .eq. '4ubpa').or.  
     &( afil(1)(1:4) .eq. '4pti' .and. afil(ii)(1:5) .eq. '1g6xa').or.  
     &( afil(1)(1:6).eq.'1hdd-C' .and. afil(ii)(1:5) .eq. '2hdda').or.  
     &( afil(1)(1:6).eq.'1bg8-A' .and. afil(ii)(1:5) .eq. '1dj8a').or.  
     &( afil(1)(1:4) .eq. '1bl0' .and. afil(ii)(1:5) .eq. '1bl0a').or.  
     &( afil(1)(1:4) .eq. '1jwe' .and. afil(ii)(1:5) .eq. '1b79a').or.  
     &( afil(1)(1:4) .eq. 'smd3' .and. afil(ii)(1:5) .eq. '1d3ba').or.  
     &( afil(1)(1:4) .eq. '1beo' .and. afil(ii)(1:5) .eq. '1lria').or.  
     &( afil(1)(1:4) .eq. '1nkl' .and. afil(ii)(1:5) .eq. '1l9la').or.  
     &( afil(1)(1:4) .eq. '1dtk' .and. afil(ii)(1:5) .eq. '1g6xa').or.  
     &( afil(1)(1:4) .eq. '2ovo' .and. afil(ii)(1:5) .eq. '1r0ri').or.  
     &( afil(1)(1:4) .eq. '1bbh' .and. afil(ii)(1:5) .eq. '1e85a').or.  
     &( afil(1)(1:4) .eq. '1bbh' .and. afil(ii)(1:5) .eq. '1s5aa').or.  
     &( afil(1)(1:4) .eq. '1c2r' .and. afil(ii)(1:5) .eq. '2bh4x').or.  
     &( afil(1)(1:4) .eq. '1cau' .and. afil(ii)(1:5) .eq. '1dgwx').or.  
     &( afil(1)(1:4) .eq. '1cau' .and. afil(ii)(1:5) .eq. '1dgwy').or.  
     &( afil(1)(1:4) .eq. '1cew' .and. afil(ii)(1:5) .eq. '1cewi').or.  
     &( afil(1)(1:4) .eq. '1dxt' .and. afil(ii)(1:5) .eq. '1irdb').or.  
     &( afil(1)(1:4) .eq. '1eaf' .and. afil(ii)(1:5) .eq. '1scza').or.  
     &( afil(1)(1:4) .eq. '1gky' .and. afil(ii)(1:5) .eq. '2an9a').or.  
     &( afil(1)(1:4) .eq. '1lga' .and. afil(ii)(1:5) .eq. '2boqa').or.  
     &( afil(1)(1:4) .eq. '2afn' .and. afil(ii)(1:5) .eq. '2bw4a').or.  
     &( afil(1)(1:4) .eq. '2fbj' .and. afil(ii)(1:5) .eq. '1mexl').or.  
     &( afil(1)(1:4) .eq. '2mta' .and. afil(ii)(1:5) .eq. '2c8sa').or.  
     &( afil(1)(1:4) .eq. '4sbv' .and. afil(ii)(1:5) .eq. '1c8na').or.  
     &( afil(1)(1:5).eq.'1b4bA' .and. afil(ii)(1:5) .eq. '1xxaa') .or.  
     &( afil(1)(1:5).eq. '1b72A' .and. afil(ii)(1:5) .eq. '2hdda').or.  
     &( afil(1)(1:5).eq. '1cewl' .and. afil(ii)(1:5) .eq. '1cewi').or.  
     &( afil(1)(1:5).eq. '1cqkA' .and. afil(ii)(1:5) .eq. '1oqoa').or.  
     &( afil(1)(1:5).eq. '1cy5A' .and. afil(ii)(1:5) .eq. '1cy5a').or.  
     &( afil(1)(1:5).eq. '1di2A' .and. afil(ii)(1:5) .eq. '1di2a').or.  
     &( afil(1)(1:5).eq. '1dtjA' .and. afil(ii)(1:5) .eq. '1zzka').or.  
     &( afil(1)(1:5).eq. '1g1cA' .and. afil(ii)(1:5) .eq. '1u2ha').or.  
     &( afil(1)(1:5).eq. '1gnuA' .and. afil(ii)(1:5) .eq. '1gnua').or.  
     &( afil(1)(1:5).eq. '1gpt_' .and. afil(ii)(1:5) .eq. '4sgbi').or.  
     &( afil(1)(1:5).eq. '1gyvA' .and. afil(ii)(1:5) .eq. '2a7ba').or.  
     &( afil(1)(1:5).eq. '1hbkA' .and. afil(ii)(1:5) .eq. '1hbka').or.  
     &( afil(1)(1:5).eq. '1jnuA' .and. afil(ii)(1:5) .eq. '1n91a').or.  
     &( afil(1)(1:5).eq. '1mn8A' .and. afil(ii)(1:5) .eq. '1mn8a').or.  
     &( afil(1)(1:5).eq. '1no5A' .and. afil(ii)(1:5) .eq. '1no5a').or.  
     &( afil(1)(1:5).eq. '1npsA' .and. afil(ii)(1:5) .eq. '1npsa').or.  
     &( afil(1)(1:6).eq. '1o2fB_' .and. afil(ii)(1:5).eq. '1svfa').or.  
     &( afil(1)(1:5).eq. '1orgA' .and. afil(ii)(1:5) .eq. '1ow4a')   
     &    )
     &    ) then
             goto 3474     
          endif

c          print*, 'Normal'

          filename=afil(ii)

         
          if(ii.eq.1) then

      open(11,file=fbase(1:ips2)//'/2640total/'//file0(1:ips3),
     &  status='old')
          do j=1,100
             if(filename(j:j).eq.'.') idot=j-1
          enddo

      open(12,file=fbase(1:ips2)//'/2640dep/'//file0(1:ips3)
     & //'_res.dep',status='old')


          else
      open(11,file= fbase(1:ips2)//'/2640total/'//filename(1:5),
     &  status='old')
      open(12,file= fbase(1:ips2)//'/2640dep/'//filename(1:5)//
     &'_res.dep', status='old')
          endif


       ime=-1 
       id=0
       id2=0
        numid = 0
        numhl=0
        numst=0
        ires=0
        ares0='     '
             
15      read(11,'(a80)',err=10,end=180) whole1(id)
       if(whole1(id)(1:3).eq.'TER') goto 10

       if(whole1(id)(1:4).ne.'ATOM') goto 15

       if( (
     &(whole1(id)(14:15) .eq. 'N ').or.
     &(whole1(id)(14:15) .eq. 'C ').or.
     &(whole1(id)(13:13) .eq. 'H') .or.
     &(whole1(id)(14:14) .eq. 'H') .or.
     &(whole1(id)(14:15) .eq. 'O ')
     &         ) ) goto 15



      read(whole1(id)(18:21),'(a4)') resname 
      READ (whole1(id)(31:54),'(3F8.3)') xptmp, yptmp,zptmp
      READ (whole1(id)(23:27),'(a5)') ares1

!       print*, 'resname',resname, xptmp,yptmp,zptmp, ares1
       if( whole1(id)(14:16) .eq. 'CA ') then

c         print*, 'point 0'
          idtemp = id
          id2=id2+1
          ires=ires+1
!          READ (whole1(id)(31:54),'(3F8.3)') xp(id2),yp(id2),zp(id2)
           resnuml(ires) =id2

           do j1=1,20
            if(resn(j1).eq.resname(1:3)) then
              restyp(id2)=j1
!              print*, resname, resn(j1), j1, id2, restyp(id2), 
!     &                resmap(restyp(id2))
            endif
           enddo
           idtemp=idtemp+1
           xxp(idtemp) = xptmp
           yyp(idtemp) = yptmp
           zzp(idtemp) = zptmp

           xp(idtemp) = xxp(idtemp)
           yp(idtemp) = yyp(idtemp)
           zp(idtemp) = zzp(idtemp)
!         print*, resname, xp(id2), yp(id2),zp(id2)
!         print*, 'id2',id2,resmap(restyp(id2))
!         print*, (resmap(restyp(l3)),l3=1,id2)
       endif

       if( (whole1(id)(18:20) .eq. 'GLY' .and. whole1(id)(14:15) .eq.
     &   'CA') .or. 
     &  (whole1(id)(18:20) .ne. 'GLY' .and. whole1(id)(14:15) .ne.
     &   'CA')  ) then

c    &  (whole1(id)(18:20) .ne. 'GLY' .and. whole1(id)(14:15) .ne.
c    &   'CA')  
!          print*, whole1(id)(14:16)
!       print*, 'id',id

        if (ares1.ne.ares0 .and. id .eq. 0 .and. numid .eq. 0) then
c         print*, 'point 1'
           ares0=ares1
           aresbefore = ares1
           numid=1
           xpmean=xptmp
           ypmean=yptmp
           zpmean=zptmp
        elseif (ares1.ne.ares0 .and. ares0 .eq. aresbefore .and. numid
     &   .ne. 0) then
           ares0=ares1

           aresbefore = ares1
           id=id+1
           xxp(id) = xpmean/numid
           yyp(id) = ypmean/numid
           zzp(id) = zpmean/numid

           xp(id) = xxp(id)
           yp(id) = yyp(id)
           zp(id) = zzp(id)

           xpmean=0.0
           ypmean=0.0
           zpmean=0.0

           numid=1
           xpmean=xptmp
           ypmean=yptmp
           zpmean=zptmp
       elseif (ares1 .eq. ares0 .and. ares1 .eq. aresbefore .and. numid
     &   .ne. 0) then
           numid = numid +1
           aresbefore = ares1
           xpmean = xpmean+xptmp
           ypmean = ypmean+yptmp
           zpmean = zpmean+zptmp
       endif

       endif

       natom=id
       nres=id2
       rnum=ires
      goto 15

180        ares0=ares1

           aresbefore = ares1
           id=id+1
           xxp(id) = xpmean/numid
           yyp(id) = ypmean/numid
           zzp(id) = zpmean/numid

           xp(id) = xxp(id)
           yp(id) = yyp(id)
           zp(id) = zzp(id)

       natom=id
       nres=id2
       rnum=ires
 10   continue
       close(11)

           do i=1,rnum
              read(12,*,err=20,end=20) j,tmp
              resdep(i)=exp(-tmp*dd)
c              print*, i,'  ', tmp, '  ',resdep(i)
           enddo
 20        continue
           close(12)

c        print*, 'id',id,'id2',id2,'rnum', rnum,'natom',natom

        if(ii.eq.1) then
         rnumme=rnum
         natome=natom
         do i=1,nres
         restyp1(i)=restyp(i)
         resnumt(i) = resnuml(i)
         xp0(i)=xp(i)
         yp0(i)=yp(i)
         zp0(i)=zp(i)
!         whole10(i)=whole1(i)
         enddo

         do i=1,natom
         xxp0(i)=xxp(i)
         yyp0(i)=yyp(i)
         zzp0(i)=zzp(i)
         mydep(i)=resdep(i)
!         print*, 'xxp0', xxp0(i)
         enddo

         write(*,'(70a)') (resmap(restyp1(i)),i=1,rnumme)
        else 


        do isf=1,rnumme
           if ( (isf.eq.1) .or. (isf.eq.2) .or. (isf.eq.3) .or. 
     &          (isf.eq.4) ) then
             do j=1,frgwind
              lr0=j
              nn=3*(j-1)+1
              xx(nn)=xp0(lr0)
              xx(nn+1)=yp0(lr0)
              xx(nn+2)=zp0(lr0)
              tmp_resnumt(j) = resnumt(j)   
              xxresnt(j) = resn(restyp1(j))
             enddo
           elseif ( (isf .eq. rnumme) .or. (isf .eq. (rnumme-1)) .or. 
     &           (isf .eq. (rnumme-2)) .or. (isf .eq. (rnumme-3)) ) then
             do j=rnumme-frgwind+1, rnumme
              lr0=j
              nn=3*(j-rnumme+frgwind-1)+1
              xx(nn)=xp0(lr0)
              xx(nn+1)=yp0(lr0)
              xx(nn+2)=zp0(lr0)
              tmp_resnumt(j-rnumme+frgwind) = resnumt(j)
              xxresnt(j-rnumme+frgwind) = resn(restyp1(j))
             enddo
           elseif ( (isf .ge. 5) .and. (isf .le. (rnumme-4)) ) then
             do j=-4,4
             lr0=j+isf
             nn=3*(j+4)+1
             xx(nn)=xp0(lr0)
             xx(nn+1)=yp0(lr0)
             xx(nn+2)=zp0(lr0)
             tmp_resnumt(j+4+1) = resnumt(j+isf)
             xxresnt(j+4+1) = resn(restyp1(j+isf))
             enddo
          endif 
 
!       if(ii .eq. 2 .and. isf .eq. 50 ) then
!        print*, xx(1),xx(2),xx(3),xx(4)
!       endif

!         write(*,*) 'centert', centert(1,isf),centert(2,isf),
!     &               centert(3,isf)

!         if(ii .eq. 7 .and. isf .eq. 50) then
!            print*, 'numt ', numt(isf)
!            print*, 'sphxt sphyt sphzt '
!          do k=1,numt(isf)
!           print*, sphxt(k,isf),sphyt(k,isf),sphzt(k,isf)
!          enddo
!         endif

!         if(ii .eq. 7 .and. isf .eq. 50) then
!            print*, 'center x y z '
!           print*, centert(1,isf),centert(2,isf),centert(3,isf)
!         endif

        do i=1,rnum
         disp=0.0
         dispa=0.0

c         print*, 'point 1.1'
           if ( (i.eq.1) .or. (i.eq.2) .or. (i.eq.3) .or. 
     &          (i.eq.4) ) then
c         print*, 'point 1.13'
             do j=1,frgwind
              lr0=j
              nn=3*(j-1)+1
              yy(nn)=xp(lr0)
              yy(nn+1)=yp(lr0)
              yy(nn+2)=zp(lr0)
              tmp_resnuml(j) = resnuml(j)
              yyresnl(j) = resn(restyp(j))
             enddo
           elseif ( (i .eq. rnum) .or. (i .eq. (rnum-1)) .or. 
     &           (i .eq. (rnum-2) ) .or. (i.eq. (rnum-3)) ) then
c         print*, 'point 1.16'
             do j=rnum-frgwind+1, rnum
              lr0=j
              nn=3*(j-rnum+fragwind-1)+1
              yy(nn)=xp(lr0)
              yy(nn+1)=yp(lr0)
              yy(nn+2)=zp(lr0)
              tmp_resnuml(j-rnum+frgwind) = resnuml(j)
              yyresnl(j-rnum+frgwind) = resn(restyp(j))
             enddo
           elseif ( (i .ge. 5) .and. (i .le. (rnum-4)) ) then
c         print*, 'point 1.19'
             do j=-4,4,1
             lr0=j+i
             nn=3*(j+4)+1
             yy(nn)=xp(lr0)
             yy(nn+1)=yp(lr0)
             yy(nn+2)=zp(lr0)
             tmp_resnuml(j+4+1) = resnuml(j+i)
             yyresnl(j+4+1) = resn(restyp(j+i))
             enddo
          endif 

c         print*, 'point 1.2'
!         if( ii .eq. 7 .and. (isf .eq. 50) .and. (i .eq. 236)) then
!         print*, 'before'
!         print*, 'isf', isf, 'i', i 
!         print*, 'whole1_tsample', (whole1_tsample(k),k=1,25)
!         print*, 'whole1', (whole1(i+k-1),k=1,25)
!         print*, 'whole1_lsam', (whole1_lsam(k),k=1,25)
!         endif
 

!         if(ii .eq. 7 .and.  (isf .eq. 50) .and.
!     &       i .eq. 236) then
!            print*, 'isf',isf,'i',i
!            print*, 'xpl ypl zpl '
!           do k=1,numl
!            print*, xpl(k),ypl(k),zpl(k)
!          enddo
!         endif

!         if(ii .eq. 7 .and.  (isf .eq. 50) .and.
!     &       i .eq. 236) then
!            print*, 'center x y z '
!            print*, centerl(1,236),centerl(2,236),centerl(3,236)
!         endif


!          write(*,*) 'numl',numl   

           tmpt(1) = 0
           tmpt(2) = 0
           tmpt(3) = 0
           tmpl(1) = 0
           tmpl(2) = 0
           tmpl(3) = 0

             do j=1,frgwind
              nn=3*(j-1)+1
              txx(nn)=xx(nn)-xp0(isf)
              txx(nn+1)=xx(nn+1)-yp0(isf)
              txx(nn+2)=xx(nn+2)-zp0(isf)
              tyy(nn)=yy(nn)-xp(i)
              tyy(nn+1)=yy(nn+1)-yp(i)
              tyy(nn+2)=yy(nn+2)-zp(i)
             enddo
c         print*, 'point 1.3'
!      if(ii .eq. 2 .and. isf .eq. 50 .and. i .eq. 50) then
!       print*, 'xp0',xp0(isf),'yp0',yp0(isf)
!       print*, 'xp',xp(isf),'yp',yp(isf)
!       print*, xx(1),xx(2),xx(3),xx(4)
!       print*, yy(1),yy(2),yy(3),yy(4)
!       print*, txx(1),txx(2),txx(3),txx(4)
!       print*, tyy(1),tyy(2),tyy(3),tyy(4)
!       endif
!       print*, 'point 1'
!       if(ii .eq. 2 .and. isf .eq. 50 .and. i .lt. 5) then 
!        print*, 'before fitsq i',i,xx(1),xx(2),xx(3),xx(4)
!       endif
        call fitsq2(disp,txx,tyy,tmpt,tmpl,frgwind,tt,bb,convg)

!       if(ii .eq. 2 .and. isf .eq. 50 .and. i .lt. 5) then 
!        print*, 'after fitsq i',i,xx(1),xx(2),xx(3),xx(4)
!       endif
!        print*,'disp', disp

c         print*, 'point 1.4'

        do k=1, 37
           resnl(k) = '   ' 
           xpl(k) = 0.0
           ypl(k) = 0.0
           zpl(k) = 0.0
        enddo

        do k=1, 37
           resnt(k) = '   '
        enddo

        do k=1, 37
           dist(k) = 0.0
           disl(k) = 0.0
           map(k) = .FALSE.
           sel_resnuml(k) = 0
        enddo
        cent_tmp_l = '   '



!        print*, 'before tmp_resnuml'
!        print*, (tmp_resnuml(k),k=1,9)
c       print*, 'point 2'
!       if(ii .eq. 2 .and. isf .eq. 50 .and. i .lt. 5) then 
!        print*, 'before select i',i,xx(1),xx(2),xx(3),xx(4)
!       endif
      call select_spheretl(xxp0,yyp0,zzp0,xxp,yyp,zzp,
     & xp0,yp0,zp0,xp,yp,zp,  
     & centert,xx,
     &  centerl,yy,sphxt,sphyt,sphzt,tmpt,xpl,ypl,zpl,tmpl,dist,disl,
     &ii,isf,numt,natome,rnumme,restyp1, i, numl,natom, restyp,frgwind,
     & i_count, tmp_resnumt,tmp_resnuml,resnumt, resnuml, sel_resnuml,
     & resnum_thold,  
     & resconta_thold,
     &  tmpl_name, pdbt, resnt,xxresnt,
!     &  whole10,whole1_tsample,
!     &    whole1_thold, 
     &   dbf_name, pdbl, resnl,yyresnl,
!     &    whole1,whole1_lsam, whole1_lsample,
     &    cent_thold,cent_tmp_l,map, radius)

!       if(ii .eq. 2 .and. isf .eq. 50 .and. i .lt. 5) then 
!        print*, 'after select i',i,xx(1),xx(2),xx(3),xx(4)
!       endif

!       if(numt(isf) .gt. 100 .or. numl .gt. 100) then
!        print*, 'numt', numt(isf), 'numl', numl
!       endif

!      if(i .eq. 2 .and. (isf .eq. 1 .or. isf .eq. 2 .or.
!     &                    isf .eq. 3)) then
!        print*, 'after select_'
!        print*, 'isf', isf, (resnum_thold(k,isf),k=1,20)
!      endif

!         if( ii .eq. 7 .and. (isf .eq. 50) .and. (i .eq. 236)) then
!         print*, 'after'
!         print*, 'isf', isf, 'i', i 
!         print*, 'numt', numt(isf), 'numl', numl
!         print*, 'sphxt', (sphxt(k,isf),k=1,25)
!         print*, 'whole1_thold', (whole1_thold(k,isf),k=1,25)
!         print*, '   ' 
!         print*, 'whole1_lsample', (whole1_lsample(k),k=1,25)
!         print*, 'dist', (dist(k),k=1,25)
!         print*, 'disl', (disl(k),k=1,25)
!         print*, 'map', (map(k),k=1,25)
!         endif

!         print*, 'disp', disp


!         if(ii .eq. 2 .and. isf .eq. 50 .and. i .eq. 17) then
!            print*, 'numt ', numt(isf)
!            print*, 'sphxt sphyt sphzt '
!          do k=1,numt(isf)
!           print*, sphxt(k,isf),sphyt(k,isf),sphzt(k,isf)
!          enddo
!         endif

!         if(ii .eq. 2 .and. isf .eq. 50 .and. i .eq. 17) then
!            print*, 'center x y z '
!           print*, centert(1,isf),centert(2,isf),centert(3,isf)
!         endif

!         if( ii .eq. 2 .and. (isf .eq. 50) .and. (i .eq. 17)) then
!            print*, 'bb '
!          do k=1,3
!            print*, bb(k,1),bb(k,2),bb(k,3)
!          enddo
!         endif

!         if(ii .eq. 2 .and.  (isf .eq. 50) .and.
!     &       i .eq. 17) then
!            print*, 'isf',isf,'i',i
!            print*, 'xpl ypl zpl '
!           do k=1,numl
!            print*, xpl(k),ypl(k),zpl(k)
!          enddo
!         endif

!         if(ii .eq. 2 .and.  (isf .eq. 50) .and.
!     &       i .eq. 17) then
!            print*, 'centerl x y z '
!            print*, centerl(1,17),centerl(2,17),centerl(3,17)
!         endif


              if(ii. eq. 2 .and. isf.eq.50 .and. i.eq.17) then
                do k=1,numt(isf)
            tmp0(k,1)=bb(1,1)*(sphxt(k,isf) - centert(1,isf))+
     &               bb(1,2)*(sphyt(k,isf) - centert(2,isf))+ 
     &               bb(1,3)*(sphzt(k,isf) - centert(3,isf))
            tmp0(k,2)=bb(2,1)*(sphxt(k,isf) - centert(1,isf))+
     &               bb(2,2)*(sphyt(k,isf) - centert(2,isf))+ 
     &               bb(2,3)*(sphzt(k,isf) - centert(3,isf))
            tmp0(k,3)=bb(3,1)*(sphxt(k,isf) - centert(1,isf))+
     &               bb(3,2)*(sphyt(k,isf) - centert(2,isf))+ 
     &               bb(3,3)*(sphzt(k,isf) - centert(3,isf))
                enddo
          
                do k1=1,numl
                 tmp2(k1,1) = xpl(k1)-centerl(1,i)
                 tmp2(k1,2) = ypl(k1)-centerl(2,i)
                 tmp2(k1,3) = zpl(k1)-centerl(3,i)
                enddo

!                  print*, 'tmp0            '
!                do k=1,numt(isf)
!                  print*, tmp0(k,1), tmp0(k,2), tmp0(k,3)
!                enddo
!                  print*, '            '
!                  print*, 'tmp2           '
!                do k1=1,numl
!                  print*, tmp2(k1,1), tmp2(k1,2), tmp2(k1,3)
!                enddo

             endif


!          write(*,*) 'just before pushe1','cent_tmp_l', cent_tmp_l    

c       print*, 'point 3'
!       if(ii .eq. 2 .and. isf .eq. 50 .and. i .lt. 5) then 
!        print*, 'before pushe1 i',i,xx(1),xx(2),xx(3),xx(4)
!       endif
        call pushe1(disp,disp0,isf,i,rnumme,numl,lhold,
     &  frgwind,distl,i_count,sel_resnuml,numt,
     &  resnum_lhold,
     &  filename,resnt,resnl,
     &  cent_tmp_l,cent_lhold,poshld
     &  ,typhold, xpl,ypl,zpl,bb,tt, sphholdx,sphholdy,sphholdz
     &  ,bbhold, tthold,centerl,ctrlhold,
     &  resdep_lhold,resdep,
     &  corres1,corres2,map)
         
!       if(ii .eq. 2 .and. isf .eq. 50 .and. i .lt. 5) then 
!        print*, 'after pushe1 i',i,xx(1),xx(2),xx(3),xx(4)
!       endif
!      if(i .eq. 2 .and. (isf .eq. 1 .or. isf .eq. 2 .or.
!     &                    isf .eq. 3)) then
!        print*, 'after pushe1'
!        print*, 'isf', isf, (resnum_thold(k,isf),k=1,20)
!      endif
!         if( ii .eq. 7 .and. (isf .eq. 50) .and. (i .eq. 236)) then
!         print*, '   '
!         print*, '   '
!         print*, 'after pushe1'
!         print*, 'isf', isf, 'i', i 
!         print*, 'map', (map(k),k=1,25)
!         print*, 'distl', (distl(291,k,50),k=1,25)
!         print*, 'corres1', (corres1(291,k,50),k=1,25)
!         print*, 'corres2', (corres2(291,k,50),k=1,25)
!         endif

!           write(*,*) isf,i, disp0(i,isf)

         enddo

        enddo

        endif
        if(mod(ii,10) .eq. 0) then
       write(*,*) afil(ii),rnum,ii
       endif


 3474  continue


       enddo  ! ii

!        print*, 'after all iteration'
!        print*, 'isf', '1', (resnum_thold(k,1),k=1,20)
!        print*, 'isf', '2', (resnum_thold(k,2),k=1,20)
!        print*, 'isf', '3', (resnum_thold(k,3),k=1,20)
!        print*, 'isf', '4', (resnum_thold(k,4),k=1,20)
!        print*, 'isf', '5', (resnum_thold(k,5),k=1,20)
!        print*, 'isf', '6', (resnum_thold(k,6),k=1,20)
!        print*, 'isf', '7', (resnum_thold(k,7),k=1,20)
!        print*, 'isf', '8', (resnum_thold(k,8),k=1,20)
!        print*, 'isf', '9', (resnum_thold(k,9),k=1,20)
!        print*, 'isf', '10', (resnum_thold(k,10),k=1,20)
!        print*, 'isf', '11', (resnum_thold(k,11),k=1,20)
!        print*, 'isf', '12', (resnum_thold(k,12),k=1,20)

!         print*, 'after pushe1'
!         print*, 'resmap', (resmap(k),k=1,20)
!         print*, 'typhold', (typhold(291,k,50),k=1,25)
!         print*, 'whole1_lhold', (whole1_lhold(291,k,50),k=1,25)
!         print*, 'distl', (distl(291,k,50),k=1,25)
!         print*, 'corres1', (corres1(291,k,50),k=1,25)
!         print*, 'corres2', (corres2(291,k,50),k=1,25)
!         print*, 'before distcalc'


!         do jmax=1, 480
!          if(lhold(jmax,50) .eq. 40 ) then
!           write(*,*) 'hlc(10,50)',hlc(jmax,50),'jmax',jmax,
!     & 'i',poshld(jmax,50),'numl', lhold(jmax,50),'disp',disp0(jmax,50)
!           endif
!          enddo 


c      in this subroutine we filter the number of Ca atom in
c      the sphere shaped PDB


!          write(*,*) 'just before pushe2'    
!         print*, 'disp0'
!       do i=1,maxres
!         write(*,*) (disp0(l,i),l=1,480)
!       enddo

c       print*, 'point 4'
        call pushe2(inb,typhold,numt,lhold,numlhold,sphholdx,sphholdy,
     &  sphholdz, bbhold,tthold)

!         do jmax=1, 480
!          if(lhold(jmax,50) .eq. 56 ) then
!           write(*,*) 'hlc(10,50)',hlc(jmax,50),'jmax',jmax,
!     & 'i',poshld(jmax,50),'numl', lhold(jmax,50),'disp',disp0(jmax,50)
!           endif
!          enddo 


!         print*, 'numlhold'
!       do i=1,200
!          write(*,*) (numlhold(l,i),l=1,150)
!         print*, '                  '
!       enddo

!         print*, 'lhold'
!       do i=1,200
!         write(*,*) (lhold(l,i),l=1,480)
!       enddo

!         print*, 'lhold(numlhold(l,i),i)'
!       do i=1,200
!          do l=1,150
!          j=numlhold(l,i)
!          if(j .gt. 0) then
!             write(*,*)'i', i, 'l',l, 'j', j, numlhold(l,i), 
!     &                 lhold(numlhold(l,i),i)
!          elseif (j .eq. 0) then
!             write(*,*)'i', i, 'l',l, 'j equal zero',j, numlhold(l,i), 
!     &                 lhold(numlhold(l,i),i)
!          endif
!             print*, '                  '
!          enddo
!       enddo

!          do i=1,1990
!             if(numt(i) .gt. 0) then
!             l=1
!            do k=1,inb
!              if((lhold(k,i) <= (numt(i)+2)).and.
!     &           (lhold(k,i) >= numt(i)) .and. (l .le. 150) ) then
!              print*,'i',i, 'numt', numt(i),'k', k, 'lhold', lhold(k,i)
!                 numlhold(l,i) = k    
!              print*, 'l', l, 'numlhold', numlhold(l,i)
!                 l=l+1       
!              endif
!            enddo
!
!           endif
!          enddo


c      in this subroutine we rotate the sphere shaped PDBs and
c      calculate the distance between the two



!          write(*,*) 'just before distcalc'    
c       print*, 'point 5'
        call distcalc(inb, numt,lhold,numlhold,orderml,distl,rnumme,
     &  resnum_thold,
     &  sphxt,sphyt,sphzt,sphholdx, sphholdy, sphholdz,bbhold,
     &  tthold, centert, ctrlhold, distml,
     &  resdep_lhold,mydep,
     &  corres1, corres2)

!        print*, 'after distcalc'

c       print*, 'point 6'


c      in this routine we calculate the prequency of specific
c      the sphere shaped PDB


          misf=rnumme-frgwind+1 

          write(24,*) misf,inb,frgwind,rnumme
          do isf1=1,misf

           do i=1,inb
           disp1(i,isf1)=sqrt(disp0(i,isf1))
           write(24,*) isf1,i,disp1(i,isf1),
     &    (typhold(i,k,isf1),k=1,50),
     &    '  ',' ',poshld(i,isf1)
           enddo            
          enddo

          close(24)


c     now construct profile and get the most probable sequence
        write(91,*) rnumme,nfil
        print*, rnumme,nfil

        mindist=100000000.0

        do i=1,rnumme
            do j=1,20
               fa(j,i)=0.0
            enddo

!         do k=1,frgwind
!             l=i-k+1
!             m= 
!             if(l.gt.0.and.l.le.misf) then
!                  do n=1,inb
!                  j=typhold(n,k,l)
!                  if(j.gt.0 .and. j.le.20) then
!                      fa(j,i)=fa(j,i)+1.0
!                  endif
!                  enddo
!             endif
!         enddo

!            print*, 'distances','i', i
!            print*, cent_thold(i), (cent_lhold(orderml(n,i),i),
!     &               distml(orderml(n,i),i),n=1,20)

!            if(i .eq. 376 .or. i .eq. 364 .or. i .eq. 248 
!     &                    .or. i .eq. 27) then
!         print*,'i', i, (cent_lhold(numlhold(orderml(n,i),i),i),
!     &         n=1,20 ), (hlc(numlhold(orderml(n,i),i),i),n=1,20),
!     &         (poshld(numlhold(orderml(n,i),i),i),n=1,20)
!            endif
!            print*, '     '

c
c If I want to change the number of neighbor, I need to change the lower n
c

c             if(i .eq. 25 ) then
c              print*,'i', i, (resconta_thold(l,i),l=1,100)
c              print*, '  '
c              print*, 'end of i .eq. 25 '
c             endif

c        print*, 'before typhold of secondpro'
!              if(radius == 6.0 )then
                 inbb = 40
!              elseif (radius == 8.0)then
!                 inbb = 40
!              endif

              do l=1, 37
                if(resnum_thold(l,i) .gt. 0) then
                  m = resnum_thold(l,i) 
!                       print*,'m    ', m
c                       print*, i, l, resconta_thold(l, i)
c                     print*, 'point ll'
                  do l2=1,37
c                     print*, 'point l2l2'
                    if(resnum_thold(l2, m) == i ) then
!                       print*, i, resnum_thold(l2, m), m
!                       print*,'m   ', m,  'l2    ', l2
c                       print*, i, l, resconta_thold(l, i)
                     do n=1,inbb
c                      print*, 'point nn'
                       k=numlhold(orderml(n,m),m)
!                       print*,'k    ', k
!                       print*,'from this typhold writing'
!                       print*, (typhold(k, kk, m), kk=1,20)
!                       print*,'from this resnum_thold writing'
!                       print*, (resnum_thold(kkk,m), kkk=1,20)
c                      print*, 'point kk'
c                        do o=1,300
c                         if(distl(k,o,m) .eq. l2) then
c                           j= typhold(k, o, m) 
c                           exit
c                         endif
c                        enddo
                  
                       j= typhold(k, distl(k,l2,m), m) 
!                      print*, 'j', j
                       if (j .ge. 1 .and. j .le. 20) then
!                      print*, 'j after', j
                         temp=real(floor(real(n)/5.0)+1)
                         fa(j,i)=fa(j,i)+1.0*(1.0/temp)
!                        print*, 'fa(j,i)     ', fa(j,i)
c                        fa(j,i)=fa(j,i) + 1
                       endif

                     enddo

                    endif
                  enddo
                endif
              enddo


c             if(i .eq. 25 ) then
c              print*,'i', i, (resconta_thold(l,i),l=1,100)
c              print*, '  '
c              print*, 'end of i .eq. 25 '
c             endif



c       print*, 'point 7'
!              disttot(i)=0
!              do n=1, 40 
!                   k=numlhold(orderml(n,i),i)
!                  if(k .gt. 0) then
!                     do j=1,20 
!                       if(resn(j) .eq. cent_lhold(k,i) ) then
!                         temp=real(floor(real(n)/5.0)+1)
!                         fa(j,i)=fa(j,i)+1.0*(1.0/temp)
c                         fa(j,i)=fa(j,i)+1.0

!                       endif
!                     enddo 
!                   endif
!                   temp=real(floor(real(n)/5.0)+1)
!           disttot(i)=disttot(i)+distml(orderml(n,i),i)*(1.0/real(n))
!              enddo

!              if(i .le. (rnumme-1) .and. mindist .gt. disttot(i) ) then
!                mindist = disttot(i)
!              endif

c        print*, 'before fa2 of secondpro'

             tmp=0.0
             tmp1=-5000.0
             imax(i)=-1
            do j=1,20
              if(fa(j,i).ge.tmp1) then
                 tmp1=fa(j,i)
                 imax(i)=j
              endif
              tmp=tmp+fa(j,i)
            enddo
            tmp=1./tmp
            it=restyp1(i)

            do j=1,20
            if(tmp .eq. 0.0 .or. fa(j,i) .eq. 0.0) then
              fa2(j,i)=0.0
            else
              fa2(j,i)=fa(j,i)*tmp
            endif
            enddo

c        write(*,*)'fa(j,i)     ', (fa(j,i),j=1,20)
c        write(*,*)'fa2(j,i)     ', (fa2(j,i),j=1,20)
c        write(*,*)'imax(i)     ', imax(i)
         write(91,100) i,' ',resmap(it),' ',it,
     &   (fa(j,i)*tmp,j=1,20)
 100   format (I6, A, A, A, I3, 2X, 20F10.3) 
        enddo   !i 


c       do i=1, rnumme
c        print*, 'i equal :   ', i
c        print*, 'resnum_thold  '
c        print*, (resnum_thold(kkk,i), kkk=1,20)
c        print*, 'resconta_thold  '
c        print*, (resconta_thold(kkk,i), kkk=1,20)
c       enddo


!        do i=1,(rnumme-1)
!         temp= disttot(i) - mindist
!        enddo

          close(91)

c         do i=1,rnumme
c        write(*,*)i, 'fa(j,i)     ', (fa(j,i),j=1,20)
c        write(*,*)i, 'fa2(j,i)     ', (fa2(j,i),j=1,20)
c        write(*,*)i, 'imax(i)     ', imax(i)
c         enddo


c        print*, 'before printout of secondpro'
!        open(99,file=fbase(1:ips)//'/out_dir/'//file0(1:ips2)//'.out')
         write(*,'(a1,i5)') '>',rnumme
         write(*,'(70a)') (resmap(restyp1(i)),i=1,rnumme)
         write(*,'(70a)') (resmap(imax(i)),i=1,rnumme)
!         write(99,'(70a)') (resmap(imax(i)),i=1,rnumme)
!         close(99)
         idt = 0.0
         do i=1, rnumme
           if(resmap(restyp1(i)) .eq. resmap(imax(i)) ) then
             idt=idt+1
           endif
         enddo
         write(*,*) '      '
         idt = real(idt*100.0/rnumme)

         idt2_60=idt
         write(*,*) radius, ' idt2  ',idt2_60

      return
      end





!    #######################################################
          subroutine pushe1(disp,disp0,isf,ixf,rnumme,numl,lhold,
     &    frgwind,distl,i_count,sel_resnuml,numt,
     &    resnum_lhold, 
     &    filename,resnt,resnl,
!     &    whole1_lsample,
!     &    whole1_lhold,
     &    cent_tmp_l,cent_lhold,poshld,
     &    typhold, xpl,ypl,zpl, bb,tt, sphholdx, sphholdy, sphholdz,
     &    bbhold, tthold,centerl,ctrlhold
     &    ,resdep_lhold,resdep
     &    ,corres1,corres2,map)
          real*8 disp,disp0(480,1990)
          integer isf,ixf,rnumme,numl,lhold(480,1990),
     &            frgwind,distl(480,37,1990),
     &            i_count,sel_resnuml(37),numt(1990),
!     &            ami(20,1270),ami_lhold(480,37,20,1990),
     &            resnum_lhold(480,37,1990)
          character filename*40, resnt(37)*3, resnl(37)*3, 
c     &          hlc(480,1990)*10, 
!     &        whole1_lsample(37)*80,
!     &          whole1_lhold(480,37,1990)*80,
     &     cent_tmp_l*3, cent_lhold(480,1990)*3
          integer poshld(480, 1990), typhold(480,37, 1990)
          real*8 xpl(37), ypl(37), zpl(37), bb(3,3), tt(3),
     &    sphholdx(480,37,1990),sphholdy(480,37,1990),
     &  sphholdz(480,37,1990),bbhold(480,9,1990),tthold(480,3,1990),
     &    centerl(3,1270), ctrlhold(480,3,1990)
     &    ,resdep_lhold(480,37,1990), resdep(1990)
          logical corres1(480,37,1990), corres2(480,37,1990),
     &     map(37)
            
         

          integer inb, me, ii, j, k, jmax, j1 
          real*8 emax
          character resn(20)*3
      data resn/'CYS','MET','PHE','ILE','LEU','VAL','TRP','TYR',
     &          'ALA','GLY','THR','SER','GLN','ASN','GLU','ASP',
     &          'HIS','ARG','LYS','PRO'/

           inb = 480 
           emax=-1.d8
           me=isf


!      if((numl .le. (numt(isf)+1)) .and. (numl .ne. 0)
!     &    .and. (numl .ge. numt(isf))  ) then


        do ii=1,inb
           if(disp0(ii,me).gt.emax) then
             jmax=ii
             emax=disp0(ii,me)
           endif
        enddo

!         if((isf .eq. 50) .and. (ixf .eq. 10)) then
!            print*, 'numl ', numl
!            print*, 'xpl ypl zpl '
!          do k=1,numl
!            print*, xpl(k),ypl(k),zpl(k)
!          enddo
!         endif

!         if( (isf .eq. 50) .and. (ixf .eq. 10)) then
!            print*, 'bb '
!          do k=1,3
!            print*, bb(k,1),bb(k,2),bb(k,3)
!          enddo
!         endif

!         write(*,*) 'resnl: after', (resnl(k), k=1,30)

          if(disp.lt.emax) then 
!            if(jmax .eq. 10 .and. isf .eq. 50 .and. numl .eq. 9) then
!               write(*,*) 'filename',filename,'isf',isf,'i',ixf
!            endif
             do j=1,37
                typhold(jmax,j,me)= 0
                sphholdx(jmax,j,me)= 0.0
                sphholdy(jmax,j,me)= 0.0
                sphholdz(jmax,j,me)= 0.0
!                whole1_lhold(jmax,j,me)= ' '
                distl(jmax,j,me)= 0
                corres1(jmax,j,me)= .FALSE.
                corres2(jmax,j,me)= .FALSE.
                resnum_lhold(jmax,j,me) = 0
                resdep_lhold(jmax,j,me) = 0.0

             enddo
             ctrlhold(jmax,1,me)= 0.0
             ctrlhold(jmax,2,me)= 0.0
             ctrlhold(jmax,3,me)= 0.0
             do j=1,3
                  tthold(jmax,j,me) = 0.0
               do k=1,3
                  bbhold(jmax,(j-1)*3+k,me) = 0.0
               enddo
             enddo
!             i_cnthold(jmax,me)= 0
             lhold(jmax,me)= 0
             disp0(jmax,me)= 0.0
c             hlc(jmax,me)= ' '
             poshld(jmax,me)= 0
             cent_lhold(jmax,me)= '   '

             do j=1,numl

                  do j1=1,20
                     if(resn(j1) .eq. resnl(j)) then
                        typhold(jmax,j,me)=j1
                     endif

!                  ami_lhold(jmax,j,j1,me) = ami(j1,sel_resnuml(j))
                  enddo
                  
                sphholdx(jmax,j,me)=xpl(j)
                sphholdy(jmax,j,me)=ypl(j)
                sphholdz(jmax,j,me)=zpl(j)

                resnum_lhold(jmax,j,me)= sel_resnuml(j)
                resdep_lhold(jmax,j,me)= resdep(sel_resnuml(j))
              enddo

!            if(jmax .eq. 59 .and. me .eq. 50) then
!                   print*, 'isf',isf,'ixf',ixf
!                   print*, 'sphholdx sphholdy sphholdz '
!                do k=1,numl
!                   print*, sphholdx(59,k,50),sphholdy(59,k,50),
!     &                     sphholdz(59,k,50)
!                enddo
!                   print*, ''
!            endif


!             do j=1,frgwind
!             hold(jmax,j,me)=resdep(k)
!             enddo
             
             ctrlhold(jmax,1,me)=centerl(1,ixf)
             ctrlhold(jmax,2,me)=centerl(2,ixf)
             ctrlhold(jmax,3,me)=centerl(3,ixf)

             do j=1,3
                  tthold(jmax,j,me) = tt(j)
               do k=1,3
                  bbhold(jmax,(j-1)*3+k,me) = bb(j,k)
               enddo
             enddo

             do j =1, numl
               if(map(j) .eqv. .TRUE.) then
                distl(jmax,j,me) = j 
                corres1(jmax,j,me) = .TRUE.
                corres2(jmax,j,me) = .TRUE.
               endif
             enddo

!            if(jmax .eq. 291 .and. me .eq. 50 .and. ixf .eq. 236) then
!              print*, 'in pushe1'
!              print*, (map(k),k=1,25)
!              print*, (distl(jmax,k,me),k=1,25)
!              print*, (corres1(jmax,k,me),k=1,25)
!              print*, (corres2(jmax,k,me),k=1,25)
!            endif
!             i_cnthold(jmax,me)=i_count
             cent_lhold(jmax,me)=cent_tmp_l
             lhold(jmax,me)=numl
             disp0(jmax,me)=disp
c             hlc(jmax,me)=filename
             poshld(jmax,me)=ixf
!           write(*,*) (resnl(k), k=1,30)
!           write(*,*) 'jmax:',jmax,'me:',me,
!     & 'cent_lhold', cent_lhold(jmax, me),'cent_tmp_l', cent_tmp_l
          endif

!           do jmax=1, 480
!            if(lhold(jmax,50) .eq. 9 ) then
!            write(*,*) 'hld(10,50)',hlc(jmax,50),'isf',isf,
!     &        'i',poshld(jmax,50),'numl', lhold(jmax,50)
!            endif
!           enddo 
 
!           write(*,*) 'pushe1 well finished'
!            endif
            return
          end




      subroutine select_spheretl(xxp0,yyp0,zzp0,xxp,yyp,zzp, 
     &    xp0,yp0,zp0,xp,yp,zp, 
     &    centert,xx,
     &    centerl,yy,sphxt,sphyt,sphzt,tmpt,xpl,ypl,zpl, tmpl,
     &    dist, disl,ite,isf,numt,
     & natome, rnumme, restyp1, ii, numl,natom, restyp,frgwind,i_count,
     &tmp_resnumt,tmp_resnuml,resnumt, resnuml, sel_resnuml,
     & resnum_thold, 
     & resconta_thold,
     &    tmpl_name, pdbt, resnt,xxresnt,
!     &    whole10,whole1_tsample,
!     &    whole1_thold, 
     &     dbf_name, pdbl, resnl,yyresnl,
!     &    whole1,whole1_lsam, whole1_lsample,
     &    cent_thold, cent_tmp_l,map,radius)
      real*8 xxp0(3750),yyp0(3750),zzp0(3750),
     &  xxp(3750),yyp(3750),zzp(3750),
     &  xp0(1270),yp0(1270),zp0(1270),
     &  xp(1270),yp(1270),zp(1270),
     &  centert(3,1990), xx(30), centerl(3,1270),yy(30),
     &  sphxt(37,1990),sphyt(37,1990),sphzt(37,1990),tmpt(3),
     &  xpl(37), ypl(37), zpl(37),tmpl(3), dist(37),disl(37)
         integer ite,isf,numt(1990),natome,rnumme,restyp1(1990),ii,numl,
     &  natom, restyp(1270),frgwind,i_count,
     &tmp_resnumt(9),tmp_resnuml(9),resnumt(1990), resnuml(1270), 
     & sel_resnuml(37),
     & resnum_thold(37,1990), 
     & resconta_thold(300,1990)
         character tmpl_name*40, pdbt*40, resnt(37)*3, xxresnt(9)*3,
!     & whole10(1990)*80,whole1_tsample(37)*80,whole1_thold(37,1990)*80,
     & dbf_name*40, pdbl*40, resnl(37)*3,yyresnl(9)*3,
!     & whole1(1270)*80,whole1_lsam(27)*80, whole1_lsample(82)*80,
     & cent_thold(1990)*3, cent_tmp_l*3
         logical map(37)
         real*8 radius


         integer lr0
         integer  i,i_count1, j, k, l, l2, tempi, centr, radi
         real temp, temp1, rad , rad0
         character whole*80, tresnt*3, tresnl*3
         real xp00(10), yp00(10), zp00(10)
         real tempx, tempy, tempz
         logical todoit

          character resn(20)*3
      data resn/'CYS','MET','PHE','ILE','LEU','VAL','TRP','TYR',
     &          'ALA','GLY','THR','SER','GLN','ASN','GLU','ASP',
     &          'HIS','ARG','LYS','PRO'/

!         write(*,*) 'centr', centr, 'dbf_name', dbf_name
         rad0 = 0.00001

         do i=1, 37
          map(i)=.FALSE.
         enddo

!        print*, 'pp   0'
!        print*, 'after tmp_resnuml'

!      if(ii .eq. 2 .and. (isf .eq. 1 .or. isf .eq. 2 .or.
!     &                    isf .eq. 3)) then
!        print*, 'isf', isf, (tmp_resnumt(k),k=1,9)
!        endif

!           do k=1, natome
!            print*,'xxp0', xxp0(k)
!           enddo

        tmpt(1) = 0.0
        tmpt(2) = 0.0
        tmpt(3) = 0.0

           tmpt(1) = xp0(isf)
           tmpt(2) = yp0(isf)
           tmpt(3) = zp0(isf)
           centert(1,isf) = tmpt(1)
           centert(2,isf) = tmpt(2)
           centert(3,isf) = tmpt(3)

        tmpl(1) = 0.0
        tmpl(2) = 0.0
        tmpl(3) = 0.0

           tmpl(1) = xp(ii)
           tmpl(2) = yp(ii)
           tmpl(3) = zp(ii)
           centerl(1,ii) = tmpl(1)
           centerl(2,ii) = tmpl(2)
           centerl(3,ii) = tmpl(3)
!          write(*,*) 'center residue :', tmpl(1), tmpl(2), tmpl(3) 


!       write(*,*) 'tmpl after', tmpl_name
         cent_thold(isf) = resn(restyp1(isf))
         cent_tmp_l = resn(restyp(ii))
!         print*, 'cent_tmp_l in select', cent_tmp_l 
          i_count1 = 1
         do k=1, 9
         temp = sqrt( (xx(3*(k-1)+1)-tmpt(1))**2 + 
     &                (xx(3*(k-1)+2)-tmpt(2))**2 +
     &                (xx(3*(k-1)+3)-tmpt(3))**2 )

!         if(ite .eq. 7 .and. isf .eq. 50 .and. ii .eq. 236) then
!            print*, 'k', k, 't temp', temp
!         endif

           if(temp <= radius ) then
           sphxt(i_count1,isf) = xx(3*(k-1) + 1)
           sphyt(i_count1,isf) = xx(3*(k-1) + 2)
           sphzt(i_count1,isf) = xx(3*(k-1) + 3)
!           resnt(i_count1) = xxresnt(k)
           resnum_thold(i_count1,isf) = tmp_resnumt(k)

             do l2=1,300
              if( resconta_thold(l2, tmp_resnumt(k)) .eq. isf) then
                todoit = .FALSE.
                exit
              else 
                todoit = .TRUE.
              endif
             enddo

           if(ite .eq. 2) then
            do l=1,300

c             do l2=1,300
c              if( resconta_thold(l2, tmp_resnumt(k)) .eq. isf) then
c                todoit = .FALSE.
c              endif
c              if(todoit .eqv. .FALSE.) exit
c             enddo

             if((l .eq. 1 .and. 
     &             resconta_thold(l, tmp_resnumt(k))
     &          .eq. 0 ) .or. 
     & (l .gt. 1 .and. resconta_thold(l-1, tmp_resnumt(k)) .ne. 0 .and.
     &               resconta_thold(l, tmp_resnumt(k)) .eq. 0) ) then
c             do l2=1,300
c              if( resconta_thold(l2, tmp_resnumt(k)) .eq. isf) then
c                todoit = .FALSE.
c              endif
c              if(todoit .eqv. .FALSE.) exit
c             enddo
!                 if(todoit .eqv. .FALSE.) then
!                  print*, 'point 1 after exit','l2',l2
!                 endif
                   if(todoit .eqv. .TRUE.) then                    
                     resconta_thold(l, tmp_resnumt(k)) = isf

c                     if(isf .ge. 10 .and. isf .le. 27) then
c                     print*, 'l', l, 'tmp_resnumt', tmp_resnumt(k), 
c     &        'resconta_thold', resconta_thold(l, tmp_resnumt(k)),
c     &'resconta_thold_total', (resconta_thold(m, tmp_resnumt(k)),m=1,l) 
c                     endif

                      goto 1101
                   endif

                 endif
!              if(todoit .eqv. .TRUE.) EXIT
!                 if(todoit .eqv. .TRUE.) then
!                  print*, 'point 2 after exit'
!                 endif
            enddo
!                  print*, 'point 1 after exit','l2',l2
           endif

!          whole1_thold(i_count1,isf) = whole1_tsample(k)
 1101      continue
           dist(i_count1) = temp
           i_count1 = i_count1 + 1
           endif
         enddo




!        print*, 'pp   1'
          i_count = 1
         do k=1, 9
         temp = sqrt( (xx(3*(k-1)+1)-tmpt(1))**2 + 
     &                (xx(3*(k-1)+2)-tmpt(2))**2 +
     &                (xx(3*(k-1)+3)-tmpt(3))**2 )


         temp1 = sqrt( (yy(3*(k-1)+1)-tmpl(1))**2 + 
     &                (yy(3*(k-1)+2)-tmpl(2))**2 +
     &                (yy(3*(k-1)+3)-tmpl(3))**2 )

!         if(ite .eq. 7 .and. isf .eq. 50 .and. ii .eq. 236) then
!            print*, 'k', k, 'l temp', temp1
!         endif
           if(temp1 <= radius ) then
           xpl(i_count) = yy(3*(k-1) + 1)
           ypl(i_count) = yy(3*(k-1) + 2)
           zpl(i_count) = yy(3*(k-1) + 3)
           resnl(i_count) = yyresnl(k)
           sel_resnuml(i_count) = tmp_resnuml(k)

!           whole1_lsample(i_count) = whole1_lsam(k)
           disl(i_count) = temp1

           if(temp <= radius .and. temp1 <= radius ) then
             map(i_count) = .TRUE.
           endif
           i_count = i_count + 1
           endif 
         enddo



!        print*, 'pp   2'
!         open(2, file='new_inp2') 
       if(ite .eq. 2) then
         i=i_count1
           do k=1, natome
           temp = sqrt((xxp0(k)-tmpt(1))**2 + (yyp0(k)-tmpt(2))**2 
     &             + (zzp0(k)-tmpt(3))**2 )  
!           print*, 'xxp0',xxp0(k), 'temp', temp
           if( temp <= radius )  then
                do j=1, i_count1-1
                 if( (abs(xxp0(k)-sphxt(j,isf)) .lt. 0.0001))  then
                    goto 17
                 endif 
                enddo
!                    write(2,100) whole
                    sphxt(i,isf) = xxp0(k)
                    sphyt(i,isf) = yyp0(k)
                    sphzt(i,isf) = zzp0(k)
!                    resnt(i) = resn(restyp1(k))
                resnum_thold(i,isf) = resnumt(k)

                 dist(i) =  temp
                    i=i+1
           endif
 17       continue
          enddo

         numt(isf)  = i - 1
!         if(isf .eq. 50) then
!         print*, 'numt',numt(isf),'isf',isf
!         endif
       endif


!        print*, 'pp   2.5'

           do k=1, rnumme
           temp = sqrt((xp0(k)-tmpt(1))**2 + (yp0(k)-tmpt(2))**2
     &             + (zp0(k)-tmpt(3))**2 )

           if( temp <= radius )  then

             do l2=1,300
              if( resconta_thold(l2, resnumt(k)) .eq. isf) then
                todoit = .FALSE.
                exit
              else 
                todoit = .TRUE.
              endif
             enddo

             if(ite .eq. 2) then
              do l=1,300

              if((l .eq. 1 .and. 
     &                resconta_thold(l, resnumt(k)) .eq. 0 ) .or. 
     & (l .gt. 1 .and. resconta_thold(l-1, resnumt(k)) .ne. 0 .and.
     &                   resconta_thold(l, resnumt(k)) .eq. 0) ) then
c             do l2=1,300
c              if( resconta_thold(l2, resnumt(k)) .eq. isf) then
c                 todoit = .FALSE.
c              endif
c              if(todoit .eqv. .FALSE.) exit
c             enddo
                      if(todoit .eqv. .TRUE.) then
                       resconta_thold(l, resnumt(k)) = isf

c                     if(isf .ge. 10 .and. isf .le. 27) then
c                    print*, 'l', l, 'resnumt', resnumt(k), 
c     &        'resconta_thold', resconta_thold(l, resnumt(k)),
c     &'resconta_thold_total', (resconta_thold(m, resnumt(k)),m=1,l) 
c                     endif

                       goto 1102
                     endif

               endif

             enddo
            endif

           endif

 1102     continue
          enddo





!        print*, 'pp   3'
!         open(1, file=dbf_name) 
!         open(2, file='new_inp3') 
         i=i_count
           do k=1, natom
           temp = sqrt((xxp(k)-tmpl(1))**2 + (yyp(k)-tmpl(2))**2 
     &             + (zzp(k)-tmpl(3))**2 )  
           if( temp <= radius  )  then
               do j=1, i_count-1
                 if( (abs(xxp(k)-xpl(j)) .lt. 0.0001))  then
                      goto 18    
                 endif
               enddo
!                    write(*,*) trim(whole), 'dist: ', temp
!                    write(2,100) whole
                    xpl(i) = xxp(k)
                    ypl(i) = yyp(k) 
                    zpl(i) = zzp(k)
                    disl(i) =  temp

                    sel_resnuml(i) = resnuml(k)
                    resnl(i) = resn(restyp(k))

!                    whole1_lsample(i) = whole1(k)
                    i=i+1
           endif
18             continue 
         enddo
!20       write(2,100) 'END' 
!100      format(A)
         numl = i - 1
!         write(*,*) 'befoer numl :', numl
!         write(*,*) (resnl(k), k=1,30)
!         pdbl = 'new_inp3'
!         close(1)
!         close(2)

c      if(ii .eq. 2 .and. (isf .eq. 11 .or. isf .eq. 22 .or.
c     &                    isf .eq. 33)) then
c        print*, 'isf', isf, (resnum_thold(k,isf),k=1,20)
c        endif

c         do k=1, 9
c                     if(isf .ge. 10 .and. isf .le. 27) then
c                     print*, 'l', l, 'tmp_resnumt', tmp_resnumt(k), 
c     &        'resconta_thold', resconta_thold(l, tmp_resnumt(k)),
c     &'resconta_thold_total', (resconta_thold(m, tmp_resnumt(k)),m=1,l) 
c                     endif
c         enddo

         return 
         end



        subroutine pushe2(inb,typhold, numt,lhold,numlhold,sphholdx,
     &  sphholdy, sphholdz,bbhold,tthold)
          integer inb,typhold(480,37, 1990),numt(1990),lhold(480,1990),
     &           numlhold(150,1990)
          real*8 sphholdx(480,37,1990),sphholdy(480,37,1990),
     &    sphholdz(480,37,1990), bbhold(480,9,1990), tthold(480,3,1990)
        
                
          integer i,j,k,l 
          write(*,*) 'inside pushe2'    
          do i=1,1990
             if(numt(i) .gt. 0) then
             l=1
            do k=1,inb
!              if(i .eq. 50) then
!              print*,'i',i, 'numt', numt(i),'k', k, 'lhold', lhold(k,i)

              if(
c     &           (lhold(k,i) <= (numt(i)+1)) .and. 
     &            (lhold(k,i) .ne. 0)
     &    .and. (lhold(k,i) == numt(i)) .and. (l .le. 150) ) then
!              if((lhold(k,i) == numt(i)) .and. (lhold(k,i) .ne. 0)
!     &    .and.  (l .le. 150) ) then
                 numlhold(l,i) = k    
!              if(i .eq. 50) then
!              print*,'i',i, 'numt', numt(i),'k', k, 'lhold', lhold(k,i)
!              print*, 'l', l, 'numlhold', numlhold(l,i)
!              endif
!              print*, 'l', l, 'numlhold', numlhold(l,i)
                 l=l+1       

              endif

!              endif
            enddo

           endif
          enddo
!             do l=1, 150
!              print*, 'l', l, 'numlhold', numlhold(l,50)
!             enddo

          write(*,*) 'end of  pushe2'    
         return
         end


      subroutine distcalc(inb, numt,lhold,numlhold,orderml,distl,rnumme,
     &  resnum_thold,
     &  sphxt, sphyt,sphzt,sphholdx, sphholdy, sphholdz,
     &  bbhold,tthold, 
     &  centert, ctrlhold, distml,
     &  resdep_lhold,mydep,
     &  corres1, corres2)
        integer inb,numt(1990), lhold(480,1990),numlhold(150,1990),
     &  orderml(150,1990),distl(480,37,1990),rnumme,
     &  resnum_thold(37,1990)
!     &  i_cnthold(480,1990)
        real*8 sphxt(37,1990),sphyt(37,1990),sphzt(37,1990),
     &  sphholdx(480,37,1990),sphholdy(480,37,1990),
     &  sphholdz(480,37,1990),bbhold(480,9,1990),tthold(480,3,1990),
     &  centert(3,1990), ctrlhold(480,3,1990),distml(150,1990),
     &  resdep_lhold(480,37,1990), mydep(1990)
        logical corres1(480,37,1990), corres2(480,37,1990)
!        character hlc(480,1990)*10


        integer i,j,k,k1,k2,k3, l,temp, minhld1, minhld2
        real*8 tmp(37,3),tmp2(37,3),tmpm,tmpdep, 
     &         mintmp
        logical yes, corre1(150), corre2(150)
          
          write(*,*) 'inside  distcalc'    
!         do i=1,1990
!           do k=1,37 
!             distml(k,i)=0.0
!           enddo
!         enddo


       do i=1,rnumme
!           ordermin = 10.0
!           orderhld = 0


!         if(i .eq. 50) then
!            print*, 'numt ', numt(i)
!            print*, 'sphxt sphyt sphzt '
!          do k=1,numt(i)
!            print*, sphxt(k,i),sphyt(k,i),sphzt(k,i)
!          enddo
!         endif

!              if(i .eq. 50) then
!                   print*, 'centert x  y  z'
!                   print*, centert(1,i),centert(2,i),
!     &                   centert(3,i)
!              endif

           do l=1,150
              j=numlhold(l,i)
!              if(i .eq. 50 .and. j .gt. 0) then
!                 print*, 'jjjjj', j
!              endif


!              if(j.eq.17 .and. i .eq. 50) then
!                   print*, 'bbholdx bbholdy bbholdz '
!                   print*, 'lhold', lhold(j,50)
!                do k=1,3
!                   print*, bbhold(j,3*(k-1)+1,50),
!     &          bbhold(j,3*(k-1)+2,50), bbhold(j,3*(k-1)+3,50)
!                enddo
!            endif

!              if(j.eq.17 .and. i .eq. 50) then
!                   print*, 'l', l, 'sphholdx sphholdy sphholdx '
!                  do k=1,lhold(j,50)
!                   print*, sphholdx(j,k,i),sphholdy(j, k,i),
!     &                   sphholdz(j,k,i)
!                  enddo
!               print*, 'lhold and hlc', lhold(j,50), hlc(j,50) 
!              endif
         
!              if(j.eq.17 .and. i .eq. 50) then
!                   print*, 'centerl x  y  z'
!                   print*, ctrlhold(j,1,i),ctrlhold(j,2,i),
!     &                   ctrlhold(j,3,i)
!              endif

              do k=1,37 
                do k1=1,3
                  tmp(k,k1) = 0.0
                  tmp2(k,k1) = 0.0
                enddo
              enddo

              if(numlhold(l,i) .gt. 0) then
!              write(*,*) 'j:', j
                do k=1,numt(i)
            tmp(k,1)=bbhold(j,1,i)*(sphxt(k,i) - centert(1,i))+
     &               bbhold(j,2,i)*(sphyt(k,i) - centert(2,i))+ 
     &               bbhold(j,3,i)*(sphzt(k,i) - centert(3,i))
            tmp(k,2)=bbhold(j,4,i)*(sphxt(k,i) - centert(1,i))+
     &               bbhold(j,5,i)*(sphyt(k,i) - centert(2,i))+ 
     &               bbhold(j,6,i)*(sphzt(k,i) - centert(3,i))
            tmp(k,3)=bbhold(j,7,i)*(sphxt(k,i) - centert(1,i))+
     &               bbhold(j,8,i)*(sphyt(k,i) - centert(2,i))+ 
     &               bbhold(j,9,i)*(sphzt(k,i) - centert(3,i))
                enddo
          
                do k1=1,lhold(j,i)
                 tmp2(k1,1) = sphholdx(j,k1,i)-ctrlhold(j,1,i)
                 tmp2(k1,2) = sphholdy(j,k1,i)-ctrlhold(j,2,i)
                 tmp2(k1,3) = sphholdz(j,k1,i)-ctrlhold(j,3,i)
                enddo

!              if(j.eq.17 .and. i.eq.50) then
!                  print*, 'tmp            '
!                do k=1,numt(i)
!                  print*, tmp(k,1), tmp(k,2), tmp(k,3)
!                enddo
!                  print*, '            '
!                  print*, 'tmp2           '
!                do k1=1,lhold(j,i)
!                  print*, tmp2(k1,1), tmp2(k1,2), tmp2(k1,3)
!                enddo
!               endif

c      
c     in this step, we need to check that all k number of Ca was paired 
c     so this step need to be checked and modified again              
c      

                do k=1, numt(i)
                  mintmp = 200.0
                  minhld1 = 0
                  minhld2 = 0

!                 if(j .eq. 291 .and. i .eq. 50) then
!                    write(*,*) (corres1(j,k3,i),k3=1,numt(i) ) 
!                    write(*,*) (corres2(j,k3,i),k3=1,lhold(j,i) ) 
!                    write(*,*) (distl(j,k3,i),k3=1,numt(i) ) 
!                 endif

                 do k1=1, numt(i)
                  do k2=1, lhold(j,i)
                     if((corres1(j,k1,i) .eqv. .FALSE.) .and. 
     &                  (corres2(j,k2,i) .eqv. .FALSE.) ) then
                     tmpm = sqrt( (tmp(k1,1)-tmp2(k2,1))**2 +
     &                            (tmp(k1,2)-tmp2(k2,2))**2 +
     &                            (tmp(k1,3)-tmp2(k2,3))**2 ) 
!                    write(*,*) numt(i),'k1', k1, 
!     &              lhold(j,i),'k2', k2, 'tmpm is', tmpm
                       if(tmpm.lt.mintmp) then
                         mintmp=tmpm 
                         minhld1= k1
                         minhld2= k2
                       endif
                     endif 
                  enddo
!                  write(*,*) '         '
!                  write(*,*) 'k',k,'mintmp',mintmp,'matching',distl(k)
                enddo
!                  write(*,*) '         '
!                  write(*,*) 'k1',minhld1,'k2',minhld2,
!     &                       'mintmp',mintmp
                  if(minhld1 .ne. 0) then
                   distl(j,minhld1,i)=minhld2
                   corres1(j,minhld1,i)=.TRUE.
                   corres2(j,minhld2,i)=.TRUE.
                  endif
               enddo 

c
c         By this point, all the corresponding mapping is ended. !!
c

c       distance calculation
                tmpm=0.0
                tmpdep=0.0
                do k=1, numt(i)
                  mintmp = 10.0
                  minhld = 0

!                  tmpm = tmpm+(tmp(k,1)-tmp2(distl(j,k,i),1))**2 +
!     &                        (tmp(k,2)-tmp2(distl(j,k,i),2))**2 +
!     &                        (tmp(k,3)-tmp2(distl(j,k,i),3))**2 

                  tmpm = tmpm+((tmp(k,1)-tmp2(distl(j,k,i),1))**2 +
     &                             (tmp(k,2)-tmp2(distl(j,k,i),2))**2 +
     &     (tmp(k,3)-tmp2(distl(j,k,i),3))**2 )*1.0/real(numt(i))
                  tmpdep = tmpdep +
     &    (mydep(resnum_thold(k,i))-resdep_lhold(j,distl(j,k,i),i))**2


                enddo

!                  tmpm = sqrt(tmpm*1.0/real(numt(i)))



!                 print*, 'tmpm', tmpm

!                 if(j .eq. 291 .and. i .eq. 50) then
!                  print*, 'j eq 291 and tmpm', tmpm
!                 endif

!                 if(j .eq. 237 .and. i .eq. 50) then
!                  print*, 'j eq 237 and tmpm', tmpm
!                 endif
                 distml(l,i)=sqrt(tmpdep)
           endif

          enddo ! l
!        print*, 'end of part1 distance calculation'

!       if(i .eq. 50) then
!         write(*,'(50F10.5)') (distml(l,i),l=1,50)
!       endif

c       sorting the distance matrix
          do k1=1,150
              orderml(k1,i)=k1
          enddo
!       if(i .eq. 50) then
!        print*, 'before sorting'
!        write(*,'(30I5)') ( orderml(l,i),l=1,30)
!        write(*,'(30F10.5)') ( distml(l,i),l=1,30)
!       endif

          do k=1,150 
             corre2(k)= .FALSE.
          enddo
          temp = 0
          do k1=1,150
                mintmp = 100000.0
                minhld1 = 0
                yes = .FALSE.
           do k2 = 1,150
              if ( (distml(k1,i) .gt. 0.0 ) .and. (distml(k2,i) .gt. 
     &            0.0) .and. ( distml(k2,i) .lt. mintmp ) .and.
     &            ( corre2(k2) .eqv. .FALSE.) )  then
                 mintmp = distml(k2,i)
                 minhld1 = k2
!                 temp = orderml(k1,i)
!                 orderml(k1,i)= orderml(k2,i)
!                 orderml(k2,i)=temp
!                 orderml(k2,i)=k1
!                 orderml(k1,i)=k2
                 yes = .TRUE.
              endif
           enddo
                 if ( yes .eqv. .TRUE.) then
                    orderml(k1,i) = minhld1
                    corre2(minhld1) = .TRUE.
                 endif
!           write(*,*) (orderml(k1,i))
          enddo  

       if(i .eq. 10) then
!        write(*,*) 'after sorting         '
!        write(*,'(30I5)') ( orderml(l,i),l=1,30)
!
        write(*,'(30F10.5)') (distml(orderml(l,i),i),l=1,30)
!        write(*,*) '          '
!        write(*,*) '          '
       endif

c  now we can use orderml  j=numlhold(orderml(k,i),i)
c  check again the   numlhold process 
c  and only use the j=numlhold(orderml(k,i),i)  in typhold process  
c 
       enddo  ! i
       print*, 'end of part2 sorting'

        return
        end


      subroutine fitsq(rms,x,y,nn,t,b,non_conv)
      implicit real*8 (a-h,o-z)
c      include 'COMMON.IOUNITS'
c  x and y are the vectors of coordinates (dimensioned (3,n)) of the two
c  structures to be superimposed.  nn is 3*n, where n is the number of  
c  points.   t and b are respectively the translation vector and the    
c  rotation matrix that transforms the second set of coordinates to the 
c  frame of the first set.                                              
c  eta =  machine-specific variable                                     
                                                                        
      dimension x(3*nn),y(3*nn),t(3)                                          
      dimension b(3,3),q(3,3),r(3,3),v(3),xav(3),yav(3),e(3),c(3,3)     
      logical non_conv
      eta = z00100000                                                   
c     small=25.0*rmdcon(3)                                              
c     small=25.0*eta                                                    
c     small=25.0*10.e-10                                                
c the following is a very lenient value for 'small'                     
      small = 0.0001D0                                                  
      non_conv=.false.
      fn=nn                                                             
      do 10 i=1,3                                                       
      xav(i)=0.0D0                                                      
      yav(i)=0.0D0                                                      
      do 10 j=1,3                                                       
   10 b(j,i)=0.0D0                                                      
      nc=0                                                              
c                                                                       
      do 30 n=1,nn                                                      
      do 20 i=1,3                                                       
c      write(iout,*)'x = ',x(nc+i),'  y = ',y(nc+i)                           
      xav(i)=xav(i)+x(nc+i)/fn                                          
   20 yav(i)=yav(i)+y(nc+i)/fn                                          
   30 nc=nc+3                                                           
c                                                                       
      do i=1,3
        t(i)=yav(i)-xav(i)
      enddo

      rms=0.0d0
      do n=1,nn
        do i=1,3
          rms=rms+(y(3*(n-1)+i)-x(3*(n-1)+i)-t(i))**2
        enddo
      enddo
      rms=dabs(rms/fn)

c     write(iout,*)'xav = ',(xav(j),j=1,3)                                    
c     write(iout,*)'yav = ',(yav(j),j=1,3)                                    
c     write(iout,*)'t   = ',(t(j),j=1,3)
c     write(iout,*)'rms=',rms
      if (rms.lt.small) return
                                                                        
                                                                        
      nc=0                                                              
      rms=0.0D0                                                         
      do 50 n=1,nn                                                      
      do 40 i=1,3                                                       
      rms=rms+((x(nc+i)-xav(i))**2+(y(nc+i)-yav(i))**2)/fn              
      do 40 j=1,3                                                       
      b(j,i)=b(j,i)+(x(nc+i)-xav(i))*(y(nc+j)-yav(j))/fn                
   40 c(j,i)=b(j,i)                                                     
   50 nc=nc+3                                                           
      call sivade(b,q,r,d,non_conv)
      sn3=dsign(1.0d0,d)                                                   
      do 120 i=1,3                                                      
      do 120 j=1,3                                                      
  120 b(j,i)=-q(j,1)*r(i,1)-q(j,2)*r(i,2)-sn3*q(j,3)*r(i,3)             
      call mvvad(b,xav,yav,t)                                           
      do 130 i=1,3                                                      
      do 130 j=1,3                                                      
      rms=rms+2.0*c(j,i)*b(j,i)                                         
  130 b(j,i)=-b(j,i)                                                    
      if (dabs(rms).gt.small) go to 140                                  
*     write (6,301)                                                     
      return                                                            
  140 if (rms.gt.0.0d0) go to 150                                         
c     write (iout,303) rms                                                 
      rms=0.0d0
*     stop                                                              
c 150 write (iout,302) dsqrt(rms)                                           
  150 continue
      return                                                            
  301 format (5x,'rms deviation negligible')                            
  302 format (5x,'rms deviation ',f14.6)                                
  303 format (//,5x,'negative ms deviation - ',f14.6)                   
      end                                                               
      subroutine sivade(x,q,r,dt,non_conv)
      implicit real*8(a-h,o-z)
c  computes q,e and r such that q(t)xr = diag(e)                        
      dimension x(3,3),q(3,3),r(3,3),e(3)                               
      dimension h(3,3),p(3,3),u(3,3),d(3)                               
      logical non_conv
      eta = z00100000                                                   
      nit = 0
      small=25.0*10.e-10                                                
c     small=25.0*eta                                                    
c     small=2.0*rmdcon(3)                                               
      xnrm=0.0d0                                                          
      do 20 i=1,3                                                       
      do 10 j=1,3                                                       
      xnrm=xnrm+x(j,i)*x(j,i)                                           
      u(j,i)=0.0d0                                                        
      r(j,i)=0.0d0                                                        
   10 h(j,i)=0.0d0                                                        
      u(i,i)=1.0                                                        
   20 r(i,i)=1.0                                                        
      xnrm=dsqrt(xnrm)                                                   
      do 110 n=1,2                                                      
      xmax=0.0d0                                                          
      do 30 j=n,3                                                       
   30 if (dabs(x(j,n)).gt.xmax) xmax=dabs(x(j,n))                         
      a=0.0d0                                                             
      do 40 j=n,3                                                       
      h(j,n)=x(j,n)/xmax                                                
   40 a=a+h(j,n)*h(j,n)                                                 
      a=dsqrt(a)                                                         
      den=a*(a+dabs(h(n,n)))                                             
      d(n)=1.0/den                                                      
      h(n,n)=h(n,n)+dsign(a,h(n,n))                                      
      do 70 i=n,3                                                       
      s=0.0d0                                                             
      do 50 j=n,3                                                       
   50 s=s+h(j,n)*x(j,i)                                                 
      s=d(n)*s                                                          
      do 60 j=n,3                                                       
   60 x(j,i)=x(j,i)-s*h(j,n)                                            
   70 continue                                                          
      if (n.gt.1) go to 110                                             
      xmax=dmax1(dabs(x(1,2)),dabs(x(1,3)))                               
      h(2,3)=x(1,2)/xmax                                                
      h(3,3)=x(1,3)/xmax                                                
      a=dsqrt(h(2,3)*h(2,3)+h(3,3)*h(3,3))                               
      den=a*(a+dabs(h(2,3)))                                             
      d(3)=1.0/den                                                      
      h(2,3)=h(2,3)+sign(a,h(2,3))                                      
      do 100 i=1,3                                                      
      s=0.0d0                                                             
      do 80 j=2,3                                                       
   80 s=s+h(j,3)*x(i,j)                                                 
      s=d(3)*s                                                          
      do 90 j=2,3                                                       
   90 x(i,j)=x(i,j)-s*h(j,3)                                            
  100 continue                                                          
  110 continue                                                          
      do 130 i=1,3                                                      
      do 120 j=1,3                                                      
  120 p(j,i)=-d(1)*h(j,1)*h(i,1)                                        
  130 p(i,i)=1.0+p(i,i)                                                 
      do 140 i=2,3                                                      
      do 140 j=2,3                                                      
      u(j,i)=u(j,i)-d(2)*h(j,2)*h(i,2)                                  
  140 r(j,i)=r(j,i)-d(3)*h(j,3)*h(i,3)                                  
      call mmmul(p,u,q)                                                 
  150 np=1                                                              
      nq=1                                                              
      nit=nit+1
      if (nit.gt.10000) then
!        print '(a)','!!!! Over 10000 iterations in SIVADE!!!!!'
        non_conv=.true.
        return
      endif
      if (dabs(x(2,3)).gt.small*(dabs(x(2,2))+abs(x(3,3)))) go to 160     
      x(2,3)=0.0d0                                                        
      nq=nq+1                                                           
  160 if (dabs(x(1,2)).gt.small*(dabs(x(1,1))+dabs(x(2,2)))) go to 180     
      x(1,2)=0.0d0                                                        
      if (x(2,3).ne.0.0d0) go to 170                                      
      nq=nq+1                                                           
      go to 180                                                         
  170 np=np+1                                                           
  180 if (nq.eq.3) go to 310                                            
      npq=4-np-nq                                                       
      if (np.gt.npq) go to 230                                          
      n0=0                                                              
      do 220 n=np,npq                                                   
      nn=n+np-1                                                         
      if (dabs(x(nn,nn)).gt.small*xnrm) go to 220                        
      x(nn,nn)=0.0d0                                                      
      if (x(nn,nn+1).eq.0.0d0) go to 220                                  
      n0=n0+1                                                           
      go to (190,210,220),nn                                            
  190 do 200 j=2,3                                                      
  200 call givns(x,q,1,j)                                               
      go to 220                                                         
  210 call givns(x,q,2,3)                                               
  220 continue                                                          
      if (n0.ne.0) go to 150                                            
  230 nn=3-nq                                                           
      a=x(nn,nn)*x(nn,nn)                                               
      if (nn.gt.1) a=a+x(nn-1,nn)*x(nn-1,nn)                            
      b=x(nn+1,nn+1)*x(nn+1,nn+1)+x(nn,nn+1)*x(nn,nn+1)                 
      c=x(nn,nn)*x(nn,nn+1)                                             
      dd=0.5*(a-b)                                                      
      xn2=c*c                                                           
      rt=b-xn2/(dd+sign(dsqrt(dd*dd+xn2),dd))                            
      y=x(np,np)*x(np,np)-rt                                            
      z=x(np,np)*x(np,np+1)                                             
      do 300 n=np,nn                                                    
      if (dabs(y).lt.dabs(z)) go to 240                                   
      t=z/y                                                             
      c=1.0/dsqrt(1.0d0+t*t)                                               
      s=c*t                                                             
      go to 250                                                         
  240 t=y/z                                                             
      s=1.0/dsqrt(1.0d0+t*t)                                               
      c=s*t                                                             
  250 do 260 j=1,3                                                      
      v=x(j,n)                                                          
      w=x(j,n+1)                                                        
      x(j,n)=c*v+s*w                                                    
      x(j,n+1)=-s*v+c*w                                                 
      a=r(j,n)                                                          
      b=r(j,n+1)                                                        
      r(j,n)=c*a+s*b                                                    
  260 r(j,n+1)=-s*a+c*b                                                 
      y=x(n,n)                                                          
      z=x(n+1,n)                                                        
      if (dabs(y).lt.dabs(z)) go to 270                                   
      t=z/y                                                             
      c=1.0/dsqrt(1.0+t*t)                                               
      s=c*t                                                             
      go to 280                                                         
  270 t=y/z                                                             
      s=1.0/dsqrt(1.0+t*t)                                               
      c=s*t                                                             
  280 do 290 j=1,3                                                      
      v=x(n,j)                                                          
      w=x(n+1,j)                                                        
      a=q(j,n)                                                          
      b=q(j,n+1)                                                        
      x(n,j)=c*v+s*w                                                    
      x(n+1,j)=-s*v+c*w                                                 
      q(j,n)=c*a+s*b                                                    
  290 q(j,n+1)=-s*a+c*b                                                 
      if (n.ge.nn) go to 300                                            
      y=x(n,n+1)                                                        
      z=x(n,n+2)                                                        
  300 continue                                                          
      go to 150                                                         
  310 do 320 i=1,3                                                      
  320 e(i)=x(i,i)                                                       
      nit=0
  330 n0=0                                                              
      nit=nit+1
      if (nit.gt.10000) then
!        print '(a)','!!!! Over 10000 iterations in SIVADE!!!!!'
        non_conv=.true.
        return
      endif
      do 360 i=1,3                                                      
      if (e(i).ge.0.0d0) go to 350                                        
      e(i)=-e(i)                                                        
      do 340 j=1,3                                                      
  340 q(j,i)=-q(j,i)                                                    
  350 if (i.eq.1) go to 360                                             
      if (dabs(e(i)).lt.dabs(e(i-1))) go to 360                           
      call switch(i,1,q,r,e)                                            
      n0=n0+1                                                           
  360 continue                                                          
      if (n0.ne.0) go to 330                                            
      if (dabs(e(3)).gt.small*xnrm) go to 370                            
      e(3)=0.0d0                                                          
      if (dabs(e(2)).gt.small*xnrm) go to 370                            
      e(2)=0.0d0                                                          
  370 dt=det(q(1,1),q(1,2),q(1,3))*det(r(1,1),r(1,2),r(1,3))            
*     write (1,501) (e(i),i=1,3)                                        
      return                                                            
  501 format (/,5x,'singular values - ',3e15.5)                         
      end                                                               
      subroutine givns(a,b,m,n)                                         
      implicit real*8 (a-h,o-z)
      dimension a(3,3),b(3,3)                                           
      if (dabs(a(m,n)).lt.dabs(a(n,n))) go to 10                          
      t=a(n,n)/a(m,n)                                                   
      s=1.0/dsqrt(1.0+t*t)                                               
      c=s*t                                                             
      go to 20                                                          
   10 t=a(m,n)/a(n,n)                                                   
      c=1.0/dsqrt(1.0+t*t)                                               
      s=c*t                                                             
   20 do 30 j=1,3                                                       
      v=a(m,j)                                                          
      w=a(n,j)                                                          
      x=b(j,m)                                                          
      y=b(j,n)                                                          
      a(m,j)=c*v-s*w                                                    
      a(n,j)=s*v+c*w                                                    
      b(j,m)=c*x-s*y                                                    
   30 b(j,n)=s*x+c*y                                                    
      return                                                            
      end                                                               
      subroutine switch(n,m,u,v,d)                                      
      implicit real*8 (a-h,o-z)
      dimension u(3,3),v(3,3),d(3)                                      
      do 10 i=1,3                                                       
      tem=u(i,n)                                                        
      u(i,n)=u(i,n-1)                                                   
      u(i,n-1)=tem                                                      
      if (m.eq.0) go to 10                                              
      tem=v(i,n)                                                        
      v(i,n)=v(i,n-1)                                                   
      v(i,n-1)=tem                                                      
   10 continue                                                          
      tem=d(n)                                                          
      d(n)=d(n-1)                                                       
      d(n-1)=tem                                                        
      return                                                            
      end                                                               
      subroutine mvvad(b,xav,yav,t)                                     
      implicit real*8 (a-h,o-z)
      dimension b(3,3),xav(3),yav(3),t(3)                               
c     dimension a(3,3),b(3),c(3),d(3)                                   
c     do 10 j=1,3                                                       
c     d(j)=c(j)                                                         
c     do 10 i=1,3                                                       
c  10 d(j)=d(j)+a(j,i)*b(i)                                             
      do 10 j=1,3                                                       
      t(j)=yav(j)                                                       
      do 10 i=1,3                                                       
   10 t(j)=t(j)+b(j,i)*xav(i)                                           
      return                                                            
      end                                                               
      double precision function det (a,b,c)
      implicit real*8 (a-h,o-z)
      dimension a(3),b(3),c(3)                                          
      det=a(1)*(b(2)*c(3)-b(3)*c(2))+a(2)*(b(3)*c(1)-b(1)*c(3))         
     1  +a(3)*(b(1)*c(2)-b(2)*c(1))                                     
      return                                                            
      end                                                               
      subroutine mmmul(a,b,c)                                           
      implicit real*8 (a-h,o-z)
      dimension a(3,3),b(3,3),c(3,3)                                    
      do 10 i=1,3                                                       
      do 10 j=1,3                                                       
      c(i,j)=0.0d0                                                        
      do 10 k=1,3                                                       
   10 c(i,j)=c(i,j)+a(i,k)*b(k,j)                                       
      return                                                            
      end                                                               
      subroutine matvec(uvec,tmat,pvec,nback)                           
      implicit real*8 (a-h,o-z)
      real*8 tmat(3,3),uvec(3,nback), pvec(3,nback)                     
c                                                                       
      do 2 j=1,nback                                                    
         do 1 i=1,3                                                     
         uvec(i,j) = 0.0d0                                                
         do 1 k=1,3                                                     
    1    uvec(i,j)=uvec(i,j)+tmat(i,k)*pvec(k,j)                        
    2 continue                                                          
      return                                                            
      end                                                               



c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine xyzatm  --  single atom internal to Cartesian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "xyzatm" computes the Cartesian coordinates of a single
c     atom from its defining internal coordinate values
c
c
      subroutine xyzatm(i,ia,bond,ib,angle1,ic,angle2,chiral)
      implicit none
      real i(3),ia(3),bond,ib(3),angle1,ic(3),angle2
      integer chiral
      include 'sizes.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      real eps,rad1,rad2
      real sin1,cos1,sin2,cos2
      real cosine,sine,sine2
      real xab,yab,zab,rab
      real xba,yba,zba,rba
      real xbc,ybc,zbc,rbc
      real xac,yac,zac,rac
      real xt,yt,zt,xu,yu,zu
      real cosb,sinb,cosg,sing
      real xtmp,ztmp,a,b,c
c
c
c     convert angles to radians, and get their sines and cosines
c
      eps = 0.00000001d0
      rad1 = angle1 / radian
      rad2 = angle2 / radian
      sin1 = sin(rad1)
      cos1 = cos(rad1)
      sin2 = sin(rad2)
      cos2 = cos(rad2)
c
c     general case where the second angle is a dihedral angle
c
      if (chiral .eq. 0) then
         xab = ia(1) - ib(1)
         yab = ia(2) - ib(2)
         zab = ia(3) - ib(3)
         rab = sqrt(xab**2 + yab**2 + zab**2)
         xab = xab / rab
         yab = yab / rab
         zab = zab / rab
         xbc = ib(1) - ic(1)
         ybc = ib(2) - ic(2)
         zbc = ib(3) - ic(3)
         rbc = sqrt(xbc**2 + ybc**2 + zbc**2)
         xbc = xbc / rbc
         ybc = ybc / rbc
         zbc = zbc / rbc
         xt = zab*ybc - yab*zbc
         yt = xab*zbc - zab*xbc
         zt = yab*xbc - xab*ybc
         cosine = xab*xbc + yab*ybc + zab*zbc
         sine = sqrt(max(1.0d0-cosine**2,eps))
         if (abs(cosine) .ge. 1.0d0) then
            write (*,*) 'error'
!            write (iout,10)  i
!   10       format (/,' XYZATM  --  Undefined Dihedral',
!     &                 ' Angle at Atom',I6)
         end if
         xt = xt / sine
         yt = yt / sine
         zt = zt / sine
         xu = yt*zab - zt*yab
         yu = zt*xab - xt*zab
         zu = xt*yab - yt*xab
         i(1) = ia(1) + bond * (xu*sin1*cos2 + xt*sin1*sin2 - xab*cos1)
         i(2) = ia(2) + bond * (yu*sin1*cos2 + yt*sin1*sin2 - yab*cos1)
         i(3) = ia(3) + bond * (zu*sin1*cos2 + zt*sin1*sin2 - zab*cos1)
c
c     general case where the second angle is a bond angle
c
      else if (abs(chiral) .eq. 1) then
         xab = ib(1) - ia(1)
         yab = ib(2) - ia(2)
         zab = ib(3) - ia(3)
         rba = sqrt(xba**2 + yba**2 + zba**2)
         xba = xba / rba
         yba = yba / rba
         zba = zba / rba
         xac = ia(1) - ic(1)
         yac = ia(2) - ic(2)
         zac = ia(3) - ic(3)
         rac = sqrt(xac**2 + yac**2 + zac**2)
         xac = xac / rac
         yac = yac / rac
         zac = zac / rac
         xt = zba*yac - yba*zac
         yt = xba*zac - zba*xac
         zt = yba*xac - xba*yac
         cosine = xba*xac + yba*yac + zba*zac
         sine2 = max(1.0d0-cosine**2,eps)
         if (abs(cosine) .ge. 1.0d0) then
!         cosine=1.0
            write (*,*) 'error'
!            write (iout,20)  i
!   20       format (/,' XYZATM  --  Defining Atoms Colinear',
!     &                 ' at Atom',i6)
         end if
!         cosine=1.0
         a = (-cos2 - cosine*cos1) / sine2
         b = (cos1 + cosine*cos2) / sine2
         c = (1.0d0 + a*cos2 - b*cos1) / sine2
         if (c .gt. eps) then
            c = chiral * sqrt(c)
         else if (c .lt. -eps) then
            c = sqrt((a*xac+b*xba)**2 + (a*yac+b*yba)**2
     &                       + (a*zac+b*zba)**2)
            a = a / c
            b = b / c
            c = 0.0d0
            if (debug) then
            write (*,*) 'error'
!               write (iout,30)  ia
!   30          format (/,' XYZATM  --  Sum of Bond Angles',
!     &                    ' Too Large at Atom',i6)
            end if
         else
            c = 0.0d0
         end if
         i(1) = ia(1) + bond * (a*xac + b*xba + c*xt)
         i(2) = ia(2) + bond * (a*yac + b*yba + c*yt)
         i(3) = ia(3) + bond * (a*zac + b*zba + c*zt)
      end if
      return
      end



