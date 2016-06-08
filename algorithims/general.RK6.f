      program test

c     need to find way to keep ivarlen different for different cases

c     program 1) van der Pol (Hairer II, pp 403)
c     program 2) Pureshi and Russo 
c     program 3) Dekker 7.5.2 pp. 215 (Kaps problem   : Index 1)
c     program 4) Dekker 7.5.1 pp. 214 (Kreiss' problem: Index 2)

      parameter(is=9,ivarlen=4)
      parameter(isamp=71,jmax=81,jactual=81)

      parameter(rho=28.,sigma=10.,beta=8./3.)	!Case 5 paramters - Lorenz

      dimension aE(is,is),aI(is,is),bE(is),bI(is),cE(is),cI(is)
      dimension bEH(is),bIH(is),bD(is,4)
      dimension stageE(is),stageI(is),maxiter(is)
      dimension svpB(is,is,0:is),alpha(is,is)
      dimension al3N(is,is),al3D(is,is),bint(is)
      dimension al4N(is,is),al4D(is,is)

      dimension ustage(ivarlen,is),res(ivarlen,is)
      dimension predvec(ivarlen,is)

      dimension uvec(ivarlen),uveco(ivarlen),usum(ivarlen)
      dimension uveciter(ivarlen),uorig(ivarlen) 
      dimension uexact(ivarlen),Rnewton(ivarlen)
      dimension errvec(ivarlen),errvecT(ivarlen),tmpvec(ivarlen)
      dimension xjacinv(ivarlen,ivarlen)

      dimension cost(isamp),error1(isamp),error2(isamp),sig(isamp)
      dimension error1P(isamp),error2P(isamp)

      dimension epsave(jmax),b1save(jmax),b2save(jmax)
      dimension b1Psave(jmax),b2Psave(jmax)
      dimension epsilon(30)

      write(*,*)'what is ipred?'
      read(*,*)ipred

      do icase = 1,1                                     !   begin algorithms loop

        icount = 0        !  cost counters
        jcount = 0        !  cost counters

        do i = 1,nrk
          stageE(i) = 0.0
          stageI(i) = 0.0
          maxiter(i)= 0
        enddo
      
        call rungeadd(aE,aI,bE,bI,cE,cI,nrk,bEH,bIH,icase,bD,
     &   svpB(1,1,0),alpha,al3N,al3D,al4N,al4D,ipred)
 
        do iprob = 1,1                                 !   begin problems loop
      
c         loop over different values of stiffness epsilon 

          do jepsil = 1,jactual,1                              !  begin stiffness epsilon loop

            itmp = 11 - jmax/jactual
            ep = 1./10**((jepsil-1)/(itmp*1.))                           !  used for 81 values of ep

            write(50,*)'zone T = "ep = ',ep,'",'
            write(51,*)'zone T = "ep = ',ep,'",'
            write(52,*)'zone T = "ep = ',ep,'",'
            write(53,*)'zone T = "ep = ',ep,'",'

            do iDT = 1,isamp,1                                 !  timestep loop for vdP, Kaps, etc
c           do iDT = 1,1,1                                     ! use to determine the exact solution

            call INIT(uvec,uexact,dt,iDT,tfinal,ep,nvecLen,iprob)     !  initialize problem information

            dto = dt
            t = 0.

            do ival = 1,nvecLen
              errvecT(ival) = 0.0
            enddo

            do i = 1,nrk                                       !  initialize stage value preditor
              do ival = 1,nvecLen                              !  with trivial guess for first stage
                predvec(ival,i) = uvec(ival)
              enddo
            enddo


            do ktime = 1,100000000                      ! advance solution (time advancement loop )
      
              if(t+dt.gt.tfinal)dt = tfinal-t + 1.0e-11

              do ival = 1,nvecLen
                uveco(ival) = uvec(ival)
              enddo

              jcount = jcount + (nrk-1)                         !  keep track of total RK stages 

              do L = 1,nrk                                      !  begin RK loop

                ttI = t + cI(L)*dt
                
                call RHS(uvec,res(1,L),dt,ep,ttI,iprob)

                do ival = 1,nvecLen                             !  write the solution into a storage register
                  ustage(ival,L) = uvec(ival)
                enddo

              if(L.ne.nrk)then

                ttI = t + cI(L+1)*dt

                do ival = 1,nvecLen
                  usum(ival) = uveco(ival)
                enddo
                do LL = 1,L 
                  do ival = 1,nvecLen
                    usum(ival) = usum(ival) + aI(L+1,LL)*res(ival,LL)
                  enddo
                enddo

                do ival = 1,nvecLen
                  if(ipred.eq.2)predvec(ival,L+1)=uvec(ival)    ! previous guess as starter
                  uvec(ival)  = predvec(ival,L+1)               ! put predict into start guess
                  uorig(ival) = uvec(ival)                      ! put predict into storage for testing
                enddo

                if(L.gt.1.and.L.lt.nrk)then
                  do ival = 1,nvecLen
                    uvec(ival)  = uveco(ival)               ! put predict into start guess
                    do jpre = 2,L
                       Z = ustage(ival,jpre)-uveco(ival)
                       uvec(ival) = uvec(ival) + alpha(L+1,jpre)*Z
                       uorig(ival) = uvec(ival)               ! put predict into storage for testing
                    enddo
                  enddo
                endif

c U^{(n+1,4)} = al4_{1}*U^{(n  ,4)} + al4_{2}*U^{(n  ,5)} +
c             + al4_{3}*U^{(n  ,6)} + al4_{4}*U^{(n+1,2)} +
c             + al4_{5}*U^{(n+1,3)}


c al4_{i} = \sum_{j=1}^{2*(order)} al4N(i,j)*r^{j} /
c           \sum_{j=1}^{2*(order)} al4D(i,j)*r^{j}

                if(L.eq.2.and.ktime.ne.1)then
                  the1 = 1.
                  the2 = the1*the
                  the3 = the2*the
                  the4 = the3*the
                  the5 = the4*the
                  the6 = the5*the
                  the7 = the6*the
                  the8 = the7*the
                  do LL = 1,5
                    xnum = al3N(LL,1)*the1
     &                   + al3N(LL,2)*the2
     &                   + al3N(LL,3)*the3
     &                   + al3N(LL,4)*the4
     &                   + al3N(LL,5)*the5
     &                   + al3N(LL,6)*the6
     &                   + al3N(LL,7)*the7
     &                   + al3N(LL,8)*the8
                    xden = al3D(LL,1)*the1
     &                   + al3D(LL,2)*the2
     &                   + al3D(LL,3)*the3
     &                   + al3D(LL,4)*the4
     &                   + al3D(LL,5)*the5
     &                   + al3D(LL,6)*the6
     &                   + al3D(LL,7)*the7
     &                   + al3D(LL,8)*the8
                    bint(LL) = xnum/xden
                  enddo
                    bint(1) =  9.9518675264213746 
                    bint(2) =  4.8366852488953721 
                    bint(3) =-24.163405114569394 
                    bint(4) = 14.152132944153401
                    bint(5) =  0.94399768676237158
                  do ival = 1,nvecLen
                    Z1 = ustage(ival,3)-ustage(ival,3)
                    Z2 = ustage(ival,4)-ustage(ival,3)
                    Z3 = ustage(ival,5)-ustage(ival,3)
                    Z4 = ustage(ival,6)-ustage(ival,3)
                    Z5 = ustage(ival,2)-ustage(ival,3)
                    uvec(ival) = ustage(ival,3) + bint(1)*z1
     &                                          + bint(2)*z2
     &                                          + bint(3)*z3
     &                                          + bint(4)*z4
     &                                          + bint(5)*z5
                    uorig(ival) = uvec(ival)               ! put predict into storage for testing
                  enddo

c                 do ival = 1,nvecLen
c                   uvec(ival) = bint(1)*ustage(ival,3)
c    &                         + bint(2)*ustage(ival,4)
c    &                         + bint(3)*ustage(ival,5)
c    &                         + bint(4)*ustage(ival,6)
c    &                         + bint(5)*ustage(ival,2)
c                   uvec(ival) = bint(1)*ustage(ival,4)
c    &                         + bint(2)*ustage(ival,5)
c    &                         + bint(3)*ustage(ival,6)
c    &                         + bint(4)*ustage(ival,2)
c    &                         + bint(5)*ustage(ival,3)
c                   uorig(ival) = uvec(ival)               ! put predict into storage for testing
c                 enddo
                endif

                do k = 1,20

                  icount = icount + 1

                  do ival = 1,nvecLen
                    uveciter(ival) = uvec(ival)
                  enddo

                  call RHS(uvec,res(1,L+1),dt,ep,ttI,iprob)

                  do ival = 1,nvecLen
                    Rnewton(ival) =  
     &              uvec(ival)  - aI(L+1,L+1)*res(ival,L+1) - usum(ival)
                  enddo

            call JACOB(uvec,xjacinv,dt,ep,aI(L+1,L+1),ttI,iprob,nvecLen)

                  do i = 1,nvecLen
                    do j = 1,nvecLen
                      uvec(i) = uvec(i) - xjacinv(i,j)*Rnewton(j)
                    enddo
                  enddo

                  tmp = 0.0
                  do ival = 1,nvecLen
                     tmp = tmp + abs(uvec(ival)-uveciter(ival))
                  enddo
                  if(tmp.lt.1.e-12) go to 160                 !  kick out of newton iteration

                enddo
  160           continue
 
             if((k.gt.maxiter(L+1)).and.(ktime.ne.1))maxiter(L+1)=k

              do ival = 1,nvecLen                    !  write the solution into a storage register
                ustage(ival,L+1) = uvec(ival)
              enddo

              xnorm = 0.0                                     !  assess how good initial guess was
              snorm = 0.0                                     !  assess how good initial guess was
              do ival = 1,nvecLen
                tmpD = (uvec(ival) - uorig(ival))
                xnorm = xnorm + tmpD*tmpD
                snorm = snorm + uvec(ival)*uvec(ival)
              enddo                                         
              xnorm = sqrt(xnorm/snorm)
              stageE(L+1) = stageE(L+1) + xnorm
              stageI(L+1) = stageI(L+1) + 1.*k

              elseif(L.eq.nrk) then
            
                do ival = 1,nvecLen
                  uvec(ival) = uveco(ival)
                enddo

                do ival = 1,nvecLen
                  do LL = 1,nrk 
                    uvec(ival) = uvec(ival) + bI(LL)*res(ival,LL)
                  enddo
                enddo

                if(time.lt.tfinal-1.0e-11)then               !  begin predicted error

                do ival = 1,nvecLen
                  errvec(ival) = 0.0
                  do LL = 1,nrk 
                    errvec(ival) = errvec(ival) 
     &                       + dt*( (bI(LL)-bIH(LL))*res(ival,LL) )
                  enddo
                  errvec(ival) = abs(errvec(ival))
                enddo
                rat = 1.0

                endif                                        !  end predicted error

c              predict new values stage values :  about 100 different kinds

c                                           !  note that ipred=2 is accomplished elsewhere
                
               if(ipred.eq.1)then
c              begin with dense output
                do K=2,nrk
                  do ival = 1,nvecLen
                    predvec(ival,K) = uvec(ival)
                  enddo
                enddo

                elseif(ipred.eq.3)then

                do K=2,nrk
                  do ival = 1,nvecLen
                    predvec(ival,K) = uveco(ival)
                    the = (1.+cI(K)*rat)
                    do LL = 1,nrk
                      bb = bD(LL,1)*the
     &                   + bD(LL,2)*the*the
     &                   + bD(LL,3)*the*the*the
     &                   + bD(LL,4)*the*the*the*the
                    predvec(ival,K) = predvec(ival,K) + bb*res(ival,LL)
                    enddo
                  enddo
                enddo

                elseif(ipred.eq.4.or.ipred.eq.5)then

c              stage value predictors
c  U^(n+1,i+1) =                  \sum_{k=0}^{order} (etah_{i k}*r^(k)) * U^{n-1} +
c                \sum_{j=1}^{s-1} \sum_{k=0}^{order} ( BBh_{ijk}*r^(k)) * U^{n,j+1}

                do ival = 1,nvecLen
                  do inew=2,nrk
                    predvec(ival,inew) = 0.0
                    do j = 1,nrk
                      the = 1.
                      bb = svpB(inew,j,0)
     &                   + svpB(inew,j,1)*the
     &                   + svpB(inew,j,2)*the*the
     &                   + svpB(inew,j,3)*the*the*the
     &                   + svpB(inew,j,4)*the*the*the*the
                      predvec(ival,inew) = predvec(ival,inew) 
     &                                   + bb*ustage(ival,j)
                    enddo
                  enddo
                enddo

                endif


              endif

            enddo                                            ! end RK loop

            do ival = 1,nvecLen
              errvecT(ival) = errvecT(ival) + errvec(ival)
            enddo                                         

            t = t + dt

            if(t.ge.tfinal) go to 100

          enddo                                              ! end time advancement loop

  100  continue
       
          cost(iDT) = alog10((nrk-1)/dto)                    !  nrk - 1 implicit stages

          totalerror  = 0.0
          totalerrorP = 0.0
          do ival = 1,nvecLen
            tmpvec(ival) = abs(uvec(ival)-uexact(ival))
            if(tmpvec(ival).eq.0.0)tmpvec(ival)=1.0e-15
            totalerror  = totalerror  + tmpvec(ival)**2 
            totalerrorP = totalerrorP + errvecT(ival)**2 
          enddo
          totalerror = sqrt(totalerror/nvecLen)

          error1(iDT)  = alog10(tmpvec(1))
          error2(iDT)  = alog10(tmpvec(2))
          error1P(iDT) = alog10(errvecT(1))
          error2P(iDT) = alog10(errvecT(2))
  
          write(50,50)cost(iDT),error1(iDT)
          write(51,50)cost(iDT),error2(iDT)
          write(52,50)cost(iDT),error1P(iDT)
          write(53,50)cost(iDT),error2P(iDT)


         enddo                                             !  end  loop over different dt

           jsamp = 41 

c          if(icase.eq.1)then                        !  make sure that data has not hit machine precision
c            jsamp = 51      
c          elseif(icase.eq.2)then
c            jsamp = 51      
c          elseif(icase.eq.3)then
c            jsamp = 51      
c          elseif(icase.eq.4)then
c            jsamp = 51      
c          endif
  
           do i=1,isamp
             sig(i) = 0.0
           enddo
 
           call fit(cost,error1,jsamp,sig,0,a1,b1,siga1,sigb1,chi2,q)
           call fit(cost,error2,jsamp,sig,0,a2,b2,siga2,sigb2,chi2,q)
           call fit(cost,error1P,jsamp,sig,0,a3,b3,siga3,sigb3,chi2,q)
           call fit(cost,error2P,jsamp,sig,0,a4,b4,siga3,sigb3,chi2,q)
           write(*,60)ep,a1,b1,a2,b2,a3,b3,a4,b4

           epsave(jepsil) = log10(ep)
           b1save(jepsil) = -b1
           b2save(jepsil) = -b2
           b1Psave(jepsil) = -b3
           b2Psave(jepsil) = -b4

         enddo                                              !  end stiffness epsilon loop

         write(35,*)'zone T = "Diff Var: Implicit",'
         do j=1,jactual
           write(35,50)epsave(j),b1save(j)
         enddo
         write(35,*)'zone T = "Alge Var: Implicit",'
         do j=1,jactual
           write(35,50)epsave(j),b2save(j)
         enddo
         write(35,*)'zone T = "Diff Var: Predicted",'
         do j=1,jactual
           write(35,50)epsave(j),b1Psave(j)
         enddo
         write(35,*)'zone T = "Alge Var: Predicted",'
         do j=1,jactual
           write(35,50)epsave(j),b2Psave(j)
         enddo

        enddo                                            !  end problems loop

       enddo                                             !  end algorithms loop

   60  format( e12.5,1x,12(f8.3,1x))

   50  format( 10(e12.5,1x))

       tmp = 1.0*icount/jcount
       write(*,*)'average iterations per step',tmp

       do istage = 2,nrk
         stageE(istage) = (nrk-1)*stageE(istage)/jcount
         stageI(istage) = (nrk-1)*stageI(istage)/jcount
       enddo
       write(*,*)'error of initial guess is '
       do istage = 2,nrk
         write(*,*) istage,stageE(istage),maxiter(istage),stageI(istage)
       enddo

       stop
       end

      subroutine INIT(uvec,uexact,dt,iDT,tfinal,ep,nvecLen,iprob)
      parameter(is=9,ivarlen=4)
        dimension uvec(ivarlen),uexact(ivarlen)
        dimension ExactTot(81,3)

        if(iprob.eq.1)then                 !  van der Pol (Hairer II, pp 403)
c    vanderPol has 8 decades with 10 levels of epsilon per decade plus 1st pt of 9th decade
 
         open(unit=39,file='exact.vanderpol.data')

          rewind(39)
          do i=1,81
            read(39,*)ExactTot(i,1),ExactTot(i,2)
            ExactTot(i,3) = 1./10**((i-1)/(10.))                  !  used for 81 values of ep
          enddo

          do i=1,81
            diff = abs(ExactTot(i,3) - ep)
            if(diff.le.1.0e-10)then
              uexact(1) = ExactTot(i,1)
              uexact(2) = ExactTot(i,2)
              go to 1 
            endif
          enddo
   1      continue

          dt = 0.5/10**((iDT-1.)/20.)
          nvecLen = 2
          tfinal = 0.5

          uvec(1) = 2.
          uvec(2) = -0.6666654321121172

        elseif(iprob.eq.2)then                     !  Pureshi and Russo 
 
          open(unit=39,file='exact.pureschi.1.data')
c         open(unit=39,file='exact.pureschi.2.data')


          rewind(39)
          do i=1,81
            read(39,*)ExactTot(i,1),ExactTot(i,2)
            ExactTot(i,3) = 1./10**((i-1)/(10.))                  !  used for 81 values of ep
          enddo

          do i=1,81
            diff = abs(ExactTot(i,3) - ep)
            if(diff.le.1.0e-10)then
              uexact(1) = ExactTot(i,1)
              uexact(2) = ExactTot(i,2)
              go to 2
            endif
          enddo
   2      continue

          dt = 0.5/10**((iDT-1.)/20.)
          nvecLen = 2
          tfinal = 5.0

          pi = acos(-1.)

c  IC: problem 1   :  equilibrium IC
          uvec(1) = pi/2.
          uvec(2) = sin(pi/2.)
c  IC: problem 2   :  non-equilibrium IC
c         uvec(1) = pi/2.
c         uvec(2) = 1./2.

        elseif(iprob.eq.3)then                   ! Dekker 7.5.2 pp. 215 (Kaps problem   : Index 1)
 
          dt = 0.5/10**((iDT-1.)/20.)
          nvecLen = 2

          tfinal = 1.0
          tmp = exp(-tfinal)
          uexact(1) = tmp*tmp
          uexact(2) = tmp

c  IC: problem 1   :  equilibrium IC
          uvec(1) = 1.
          uvec(2) = 1.
c  IC: problem 2   :  non-equilibrium IC
c         uvec(1) = ?
c         uvec(2) = ?

        elseif(iprob.eq.4)then          ! Dekker 7.5.1 pp. 214 (Kreiss' problem: Index 2)
 
          dt = 0.5/10**((iDT-1.)/20.)
          nvecLen = 2

          tfinal = 1.0
          tmpP = exp(+tfinal   )
          tmpM = exp(-tfinal/ep)
          tmp = exp(tfinal)
          uo = 1.
          vo = 1.
          uexact(1) = uo*tmpP + vo*(tmpP - tmpM)/(1. + ep)
          uexact(2) = vo*tmpM

c  IC: problem 1   :  equilibrium IC
          uvec(1) = uo
          uvec(2) = vo
c  IC: problem 2   :  non-equilibrium IC
c         uvec(1) = ?
c         uvec(2) = ?
	elseif(iprob.eq.5)then
	
	  dt = 0.5/10**((iDT-1.)/20.)
	  nveclen = 2 !!!!!!!!!!! has something to do with jacobians, but not sure if it is relevant

	  tfinal = 1.0 !!arbitrary
	

        endif

      return
      end

c-------------------------------------------------------------
      subroutine RHS(uvec,res,dt,ep,ttI,iprob)
      parameter(is=9,ivarlen=4)
      dimension uvec(ivarlen),res(ivarlen)
         
      if    (iprob.eq.1)then
        res(1) = dt*uvec(2)
        res(2) = dt*((1-uvec(1)*uvec(1))*uvec(2) - uvec(1))/ep
      elseif(iprob.eq.2)then
        res(1) = dt*(-uvec(2))
        res(2) = dt*( uvec(1) + (sin(uvec(1)) - uvec(2))/ep)
      elseif(iprob.eq.3)then
        res(1) = dt*(-(1./ep+2.)*uvec(1) + (1./ep)*uvec(2)*uvec(2))
        res(2) = dt*(uvec(1) - uvec(2) - uvec(2)*uvec(2) )
      elseif(iprob.eq.4)then
        res(1) = dt*(uvec(1) + uvec(2)/ep)
        res(2) = dt*(        - uvec(2)/ep)
      endif

      return
      end

c--------------------------------------------------------------
      subroutine JACOB(uvec,xinv,dt,ep,akk,ttI,iprob,nvecLen)
      parameter(is=9,ivarlen=4)
      dimension uvec(ivarlen),xjac(ivarlen,ivarlen)
      dimension xinv(ivarlen,ivarlen)

   
      if(iprob.eq.1)then
        xjac(1,1) = 1.-akk*dt*(0.)
        xjac(1,2) = 0.-akk*dt*(1.)
        xjac(2,1) = 0.-akk*dt*(-2*uvec(1)*uvec(2)-1)/ep
        xjac(2,2) = 1.-akk*dt*(1-uvec(1)*uvec(1))/ep
      elseif(iprob.eq.2)then
        xjac(1,1) = 1.-akk*dt*(0.)
        xjac(1,2) = 0.-akk*dt*(0.)
        xjac(2,1) = 0.-akk*dt*(cos(uvec(1)) )/ep
        xjac(2,2) = 1.-akk*dt*(-1.)/ep
      elseif(iprob.eq.3)then
        xjac(1,1) = 1.-akk*dt*(-(1./ep))
        xjac(1,2) = 0.-akk*dt*(+1./ep)*2*uvec(2)
        xjac(2,1) = 0.-akk*dt*(0)
        xjac(2,2) = 1.-akk*dt*(0)
      elseif(iprob.eq.4)then
        xjac(1,1) = 1.-akk*dt
        xjac(1,2) = 0.-akk*dt/ep
        xjac(2,1) = 0.-akk*dt*(0)
        xjac(2,2) = 1.+akk*dt/ep
      endif

      if(nvecLen.eq.2)then
        det = (xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1))
        xinv(1,1) =  xjac(2,2)/det
        xinv(1,2) = -xjac(1,2)/det
        xinv(2,1) = -xjac(2,1)/det
        xinv(2,2) =  xjac(1,1)/det

      elseif(nvecLen.eq.3)then

        x11 = xjac(1,1)
        x12 = xjac(1,2)
        x13 = xjac(1,3)
        x21 = xjac(2,1)
        x22 = xjac(2,2)
        x23 = xjac(2,3)
        x31 = xjac(3,1)
        x32 = xjac(3,2)
        x33 = xjac(3,3)

        det = - x13*x22*x31 + x12*x23*x31 +  x13*x21*x32 
     &        - x11*x23*x32 - x12*x21*x33 +  x11*x22*x33

        detI = 1./det

        xinv(1,1) =  (- x23*x32 + x22*x33) * detI
        xinv(1,2) =  (+ x13*x32 - x12*x33) * detI
        xinv(1,3) =  (- x13*x22 + x12*x23) * detI
        xinv(2,1) =  (+ x23*x31 - x21*x33) * detI
        xinv(2,2) =  (- x13*x31 + x11*x33) * detI
        xinv(2,3) =  (+ x13*x21 - x11*x23) * detI
        xinv(3,1) =  (- x22*x31 + x21*x32) * detI
        xinv(3,2) =  (+ x12*x31 - x11*x32) * detI
        xinv(3,3) =  (- x12*x21 + x11*x22) * detI

      elseif(nvecLen.eq.4)then

        x11 = xjac(1,1)
        x12 = xjac(1,2)
        x13 = xjac(1,3)
        x14 = xjac(1,4)
        x21 = xjac(2,1)
        x22 = xjac(2,2)
        x23 = xjac(2,3)
        x24 = xjac(2,4)
        x31 = xjac(3,1)
        x32 = xjac(3,2)
        x33 = xjac(3,3)
        x34 = xjac(3,4)
        x41 = xjac(4,1)
        x42 = xjac(4,2)
        x43 = xjac(4,3)
        x44 = xjac(4,4)

        det = 
     -  (x14*x23*x32*x41 - x13*x24*x32*x41 - 
     -   x14*x22*x33*x41 + x12*x24*x33*x41 + 
     -   x13*x22*x34*x41 - x12*x23*x34*x41 - 
     -   x14*x23*x31*x42 + x13*x24*x31*x42 + 
     -   x14*x21*x33*x42 - x11*x24*x33*x42 - 
     -   x13*x21*x34*x42 + x11*x23*x34*x42 + 
     -   x14*x22*x31*x43 - x12*x24*x31*x43 - 
     -   x14*x21*x32*x43 + x11*x24*x32*x43 + 
     -   x12*x21*x34*x43 - x11*x22*x34*x43 - 
     -   x13*x22*x31*x44 + x12*x23*x31*x44 + 
     -   x13*x21*x32*x44 - x11*x23*x32*x44 - 
     -   x12*x21*x33*x44 + x11*x22*x33*x44)

        detI = 1./det

        xinv(1,1) = (
     - -(x24*x33*x42) + x23*x34*x42 + x24*x32*x43 - 
     -   x22*x34*x43  - x23*x32*x44 + x22*x33*x44  ) * detI
        xinv(1,2) = (
     -   x14*x33*x42  - x13*x34*x42 - x14*x32*x43 + 
     -   x12*x34*x43  + x13*x32*x44 - x12*x33*x44  ) * detI
        xinv(1,3) = (
     - -(x14*x23*x42) + x13*x24*x42 + x14*x22*x43 - 
     -   x12*x24*x43  - x13*x22*x44 + x12*x23*x44  ) * detI
        xinv(1,4) = (
     -   x14*x23*x32  - x13*x24*x32 - x14*x22*x33 + 
     -   x12*x24*x33  + x13*x22*x34 - x12*x23*x34  ) * detI
        xinv(2,1) = (
     -   x24*x33*x41  - x23*x34*x41 - x24*x31*x43 + 
     -   x21*x34*x43  + x23*x31*x44 - x21*x33*x44  ) * detI
        xinv(2,2) = (
     - -(x14*x33*x41) + x13*x34*x41 + x14*x31*x43 - 
     -   x11*x34*x43  - x13*x31*x44 + x11*x33*x44  ) * detI
        xinv(2,3) = (
     -   x14*x23*x41  - x13*x24*x41 - x14*x21*x43 + 
     -   x11*x24*x43  + x13*x21*x44 - x11*x23*x44  ) * detI
        xinv(2,4) = (
     - -(x14*x23*x31) + x13*x24*x31 + x14*x21*x33 - 
     -   x11*x24*x33  - x13*x21*x34 + x11*x23*x34  ) * detI
        xinv(3,1) = (
     - -(x24*x32*x41) + x22*x34*x41 + x24*x31*x42 - 
     -   x21*x34*x42  - x22*x31*x44 + x21*x32*x44  ) * detI
        xinv(3,2) = (
     -   x14*x32*x41  - x12*x34*x41 - x14*x31*x42 + 
     -   x11*x34*x42  + x12*x31*x44 - x11*x32*x44  ) * detI
        xinv(3,3) = (
     - -(x14*x22*x41) + x12*x24*x41 + x14*x21*x42 - 
     -   x11*x24*x42  - x12*x21*x44 + x11*x22*x44  ) * detI
        xinv(3,4) = (
     -   x14*x22*x31  - x12*x24*x31 - x14*x21*x32 + 
     -   x11*x24*x32  + x12*x21*x34 - x11*x22*x34  ) * detI
        xinv(4,1) = (
     -   x23*x32*x41  - x22*x33*x41 - x23*x31*x42 + 
     -   x21*x33*x42  + x22*x31*x43 - x21*x32*x43  ) * detI
        xinv(4,2) = (
     - -(x13*x32*x41) + x12*x33*x41 + x13*x31*x42 - 
     -   x11*x33*x42  - x12*x31*x43 + x11*x32*x43  ) * detI
        xinv(4,3) = (
     -   x13*x22*x41  - x12*x23*x41 - x13*x21*x42 + 
     -   x11*x23*x42  + x12*x21*x43 - x11*x22*x43  ) * detI
        xinv(4,4) = (
     - -(x13*x22*x31) + x12*x23*x31 + x13*x21*x32 - 
     -   x11*x23*x32  - x12*x21*x33 + x11*x22*x33  ) * detI

      endif

  
      return
      end

c------------------------------------------------------------------c

      subroutine rungeadd(aE,aI,bE,bI,cE,cI,ns,bEH,bIH,icase,bD,
     &                    svpB,alpha,al3N,al3D,al4N,al4D,ipred)
      parameter(is=9)

      dimension aE(is,is),bE(is),cE(is)
      dimension aI(is,is),bI(is),cI(is)
      dimension bEh(is),bIh(is)
      dimension bD(is,4),bdsum(is),alpha(is,is)
      dimension al3N(is,is),al3D(is,is)
      dimension al4N(is,is),al4D(is,is)
      dimension p(0:is,0:is),q(0:is)
      dimension svpB(is,is,0:is)

      do i=1,is
       do j=1,is
        aE(i,j) = 0.0
        aI(i,j) = 0.0
       enddo
       be(i)  = 0.
       bi(i)  = 0.
       beh(i) = 0.
       bih(i) = 0.
       ce(i)  = 0.
       ci(i)  = 0.
      enddo

       do i = 1,is
         bdsum(i) = 0.0
         do j = 1,4
           bd(i,j) = 0.0
         enddo
       enddo

      if(icase.eq.1) then     !    The scheme you're using is: ESDIRK436L2SA_ARK
    
      ns = 6 

c  The scheme coefficients, A, b, bh, and c  are given by
    
      ai(2,1) = 1.d0/4.d0
      ai(2,2) = 1.d0/4.d0
      ai(3,1) = 8611.d0/62500.d0
      ai(3,2) = -1743.d0/31250.d0
      ai(3,3) = 1.d0/4.d0
      ai(4,1) = 5012029.d0/34652500.d0
      ai(4,2) = -654441.d0/2922500.d0
      ai(4,3) = 174375.d0/388108.d0
      ai(4,4) = 1.d0/4.d0
      ai(5,1) = 15267082809.d0/155376265600.d0
      ai(5,2) = -71443401.d0/120774400.d0
      ai(5,3) = 730878875.d0/902184768.d0
      ai(5,4) = 2285395.d0/8070912.d0
      ai(5,5) = 1.d0/4.d0
      ai(6,1) = 82889.d0/524892.d0
      ai(6,3) = 15625.d0/83664.d0
      ai(6,4) = 69875.d0/102672.d0
      ai(6,5) = -2260.d0/8211.d0
      ai(6,6) = 1.d0/4.d0
    
      be(1) = 82889.d0/524892.d0
      be(2) = 0.d0/1.d0
      be(3) = 15625.d0/83664.d0
      be(4) = 69875.d0/102672.d0
      be(5) = -2260.d0/8211.d0
      be(6) = 1.d0/4.d0
    
      beh(1) = 4586570599.d0/29645900160.d0
      beh(2) = 0.d0/1.d0
      beh(3) = 178811875.d0/945068544.d0
      beh(4) = 814220225.d0/1159782912.d0
      beh(5) = -3700637.d0/11593932.d0
      beh(6) = 61727.d0/225920.d0
    
      bi(1) = 82889.d0/524892.d0
      bi(2) = 0.d0/1.d0
      bi(3) = 15625.d0/83664.d0
      bi(4) = 69875.d0/102672.d0
      bi(5) = -2260.d0/8211.d0
      bi(6) = 1.d0/4.d0
    
      bih(1) = 4586570599.d0/29645900160.d0
      bih(2) = 0.d0/1.d0
      bih(3) = 178811875.d0/945068544.d0
      bih(4) = 814220225.d0/1159782912.d0
      bih(5) = -3700637.d0/11593932.d0
      bih(6) = 61727.d0/225920.d0
    
c   The primary predictor is now given in two parts; eta and B
    
c   U^{(n+1,i+1)} = \sum_{j=0}^{order}eta_{ij}*r^{j} +  
c                 + \sum_{j=1}^{S-1}\sum_{k=0}^{order}B_{ijk}*r^{k}U^{(n,j+1)} 
    
      svpB(2,6,0) = 1.d0/1.d0
      svpB(3,6,0) = 1.d0/1.d0
      svpB(4,6,0) = 1.d0/1.d0
      svpB(5,6,0) = 1.d0/1.d0
      svpB(6,6,0) = 1.d0/1.d0
      svpB(2,2,1) = -3377.d0/126.d0
      svpB(3,2,1) = -29222231.d0/1968750.d0
      svpB(4,2,1) = -17389771859.d0/1091553750.d0
      svpB(5,2,1) = -1159773601703.d0/101965674300.d0
      svpB(6,2,1) = -147230135.d0/8267049.d0
      svpB(2,3,1) = 2714875.d0/31374.d0
      svpB(3,3,1) = 189518369.d0/3921750.d0
      svpB(4,3,1) = 116964555941.d0/2174375070.d0
      svpB(5,3,1) = 63034726731115.d0/1523367174042.d0
      svpB(6,3,1) = 42798736375.d0/686165067.d0
      svpB(2,4,1) = -570545.d0/19251.d0
      svpB(3,4,1) = -1036509259.d0/60159375.d0
      svpB(4,4,1) = -747062586151.d0/33354763875.d0
      svpB(5,4,1) = -411728784796537.d0/18694678056660.d0
      svpB(6,4,1) = -2767877890.d0/93561999.d0
      svpB(2,5,1) = -386912.d0/8211.d0
      svpB(3,5,1) = -3299806112.d0/128296875.d0
      svpB(4,5,1) = -1835611316768.d0/71132919375.d0
      svpB(5,5,1) = -1379738710208.d0/88204820025.d0
      svpB(6,5,1) = -28824358208.d0/1077472053.d0
      svpB(2,6,1) = 218.d0/5.d0
      svpB(3,6,1) = 1894838.d0/78125.d0
      svpB(4,6,1) = 1149924782.d0/43315625.d0
      svpB(5,6,1) = 239743897067.d0/12138770750.d0
      svpB(6,6,1) = 19865372.d0/656115.d0
      svpB(2,2,2) = 143443.d0/42084.d0
      svpB(3,2,2) = 988178827.d0/657562500.d0
      svpB(4,2,2) = 137848723.d0/26302500.d0
      svpB(5,2,2) = 41455027.d0/4208400.d0
      svpB(6,2,2) = 143443.d0/10521.d0
      svpB(2,3,2) = -119495125.d0/10478916.d0
      svpB(3,3,2) = -79344763.d0/15781500.d0
      svpB(4,3,2) = -918678521.d0/52394580.d0
      svpB(5,3,2) = -1381363645.d0/41915664.d0
      svpB(6,3,2) = -119495125.d0/2619729.d0
      svpB(2,4,2) = -27276895.d0/6429834.d0
      svpB(3,4,2) = -37582105931.d0/20093231250.d0
      svpB(4,4,2) = -169116749.d0/25926750.d0
      svpB(5,4,2) = -1576604531.d0/128596680.d0
      svpB(6,4,2) = -54553790.d0/3214917.d0
      svpB(2,5,2) = 511664.d0/26887.d0
      svpB(3,5,2) = 3524853296.d0/420109375.d0
      svpB(4,5,2) = 491709104.d0/16804375.d0
      svpB(5,5,2) = 36967724.d0/672175.d0
      svpB(6,5,2) = 2046656.d0/26887.d0
      svpB(2,6,2) = -38289.d0/3340.d0
      svpB(3,6,2) = -263772921.d0/52187500.d0
      svpB(4,6,2) = -36795729.d0/2087500.d0
      svpB(5,6,2) = -11065521.d0/334000.d0
      svpB(6,6,2) = -38289.d0/835.d0
      svpB(2,2,3) = 19693.d0/3006.d0
      svpB(3,2,3) = 67015279.d0/46968750.d0
      svpB(4,2,3) = 2948022407.d0/313751250.d0
      svpB(5,2,3) = 1098264805207.d0/45380980800.d0
      svpB(6,2,3) = 157544.d0/4509.d0
      svpB(2,3,3) = -153034375.d0/5239458.d0
      svpB(3,3,3) = -2007811.d0/315630.d0
      svpB(4,3,3) = -36654548645.d0/874989486.d0
      svpB(5,3,3) = -341384792666125.d0/3163961981376.d0
      svpB(6,3,3) = -1224275000.d0/7859187.d0
      svpB(2,4,3) = 15804500.d0/3214917.d0
      svpB(3,4,3) = 430261708.d0/401864625.d0
      svpB(4,4,3) = 610559444.d0/86595345.d0
      svpB(5,4,3) = 8814058860455.d0/485349589656.d0
      svpB(6,4,3) = 252872000.d0/9644751.d0
      svpB(2,5,3) = 13198400.d0/457079.d0
      svpB(3,5,3) = 1796566208.d0/285674375.d0
      svpB(4,5,3) = 79031491264.d0/1908304825.d0
      svpB(5,5,3) = 478959862.d0/4490129.d0
      svpB(6,5,3) = 211174400.d0/1371237.d0
      svpB(2,6,3) = -7071.d0/334.d0
      svpB(3,6,3) = -24062613.d0/5218750.d0
      svpB(4,6,3) = -1058521629.d0/34861250.d0
      svpB(5,6,3) = -394344713229.d0/5042331200.d0
      svpB(6,6,3) = -18856.d0/167.d0
      svpB(2,2,4) = 74623.d0/18036.d0
      svpB(3,2,4) = -22403609.d0/281812500.d0
      svpB(4,2,4) = 44562253271.d0/13177552500.d0
      svpB(5,2,4) = 926001497797.d0/68071471200.d0
      svpB(6,2,4) = 631898.d0/31563.d0
      svpB(2,3,4) = -547529875.d0/31436748.d0
      svpB(3,3,4) = -3601861.d0/47344500.d0
      svpB(4,3,4) = -460992042829.d0/26249684580.d0
      svpB(5,3,4) = -323698602127855.d0/4745942972064.d0
      svpB(6,3,4) = -777175250.d0/7859187.d0
      svpB(2,4,4) = 92010395.d0/19289502.d0
      svpB(3,4,4) = 1642447243.d0/60279693750.d0
      svpB(4,4,4) = 63148700849.d0/12989301750.d0
      svpB(5,4,4) = 17154291013966.d0/910030480605.d0
      svpB(6,4,4) = 263402720.d0/9644751.d0
      svpB(2,5,4) = 18331312.d0/1371237.d0
      svpB(3,5,4) = 5096323504.d0/21425578125.d0
      svpB(4,5,4) = 2137276086032.d0/143122861875.d0
      svpB(5,5,4) = 19239796526.d0/336759675.d0
      svpB(6,5,4) = 112833664.d0/1371237.d0
      svpB(2,6,4) = -35427.d0/3340.d0
      svpB(3,6,4) = -8040459.d0/52187500.d0
      svpB(4,6,4) = -4032640797.d0/348612500.d0
      svpB(5,6,4) = -559786648953.d0/12605828000.d0
      svpB(6,6,4) = -53486.d0/835.d0
    
c   U^{(n+1,i+1)} = \sum_{j=0}^{order}eta_{ij}*r^{j} +  
c                 + \sum_{j=1}^{S-1}\sum_{k=0}^{order}B_{ijk}*r^{k}U^{(n,j+1)} 
    
      svpB(2,1,1) = -17434754.d0/656115.d0
      svpB(3,1,1) = -151586831054.d0/10251796875.d0
      svpB(4,1,1) = -92113753426406.d0/5684006259375.d0
      svpB(5,1,1) = -85369397385244.d0/7062113095533.d0
      svpB(6,1,1) = -1593371832236.d0/86097378645.d0
      svpB(2,1,2) = 40139517.d0/8593820.d0
      svpB(3,1,2) = 3331579911.d0/1617812500.d0
      svpB(4,1,2) = 1244325027.d0/173262500.d0
      svpB(5,1,2) = 11600320413.d0/859382000.d0
      svpB(6,1,2) = 40139517.d0/2148455.d0
      svpB(2,1,3) = 146620279.d0/14609494.d0
      svpB(3,1,3) = 6011431439.d0/2750281250.d0
      svpB(4,1,3) = 708029327291.d0/49189223750.d0
      svpB(5,1,3) = 480994716412613.d0/12973918177600.d0
      svpB(6,1,3) = 1172962232.d0/21914241.d0
      svpB(2,1,4) = 2519213441.d0/438284820.d0
      svpB(3,1,4) = 3674587459.d0/82508437500.d0
      svpB(4,1,4) = 8784810383521.d0/1475676712500.d0
      svpB(5,1,4) = 384305477597752.d0/16689657452547.d0
      svpB(6,1,4) = 1217184446.d0/36523735.d0
    
c  The secondary predictor is used well into the next step
c  It uses only information from the new step 
    
c   Z_{0}^{(n+1,i)} = \sum_{j=2}^{i-1} \alpha_{ij} Z_{0}^{(n+1,j)} 
    
      alpha(3,2) = 83.d0/125.d0
      alpha(4,2) = 372.d0/175.d0
      alpha(4,3) = -775.d0/581.d0
      alpha(5,2) = 52003.d0/16272.d0
      alpha(5,3) = -65639125.d0/16206912.d0
      alpha(5,4) = 5825645.d0/6053184.d0
      alpha(6,2) = -155.d0/567.d0
      alpha(6,3) = -113125.d0/141183.d0
      alpha(6,4) = 43300.d0/173259.d0
      alpha(6,5) = 36160.d0/24633.d0
    
c   This secondary predictor uses information from the new and old step 
    
    
c   U^{(n+1,3)} = al3_{1}*U^{(n  ,3)} + al3_{2}*U^{(n  ,4)} +
c               + al3_{3}*U^{(n  ,5)} + al3_{4}*U^{(n  ,6)} +
c               + al3_{5}*U^{(n+1,2)}
    
    
c   al3_{i} = \sum_{j=1}^{2*(order)} al3N(i,j)*r^{j} / 
c             \sum_{j=1}^{2*(order)} al3D(i,j)*r^{j} 
    
    
      al3N(1,1) = -13406081696.d0/9273941175.d0
      al3N(2,1) = 27150465254752.d0/5332516175625.d0
      al3N(3,1) = -35615616927744.d0/8887526959375.d0
      al3N(4,1) = 1778562681504.d0/1932071078125.d0
      al3N(5,1) = 6889.d0/15625.d0
      al3N(1,2) = 71605868527.d0/27821823525.d0
      al3N(2,2) = -124718990615639.d0/15997548526875.d0
      al3N(3,2) = 67054977035008.d0/8887526959375.d0
      al3N(4,2) = -1406942044326.d0/1932071078125.d0
      al3N(5,2) = 857620765896.d0/1932071078125.d0
      al3N(1,3) = 11480922584.d0/3091313725.d0
      al3N(2,3) = -24748852324168.d0/1777505391875.d0
      al3N(3,3) = 27909908271744.d0/1777505391875.d0
      al3N(4,3) = -1647422231952.d0/386414215625.d0
      al3N(5,3) = -2509229108.d0/386414215625.d0
      al3N(1,4) = 8282570581.d0/3091313725.d0
      al3N(2,4) = -665004331789.d0/77282843125.d0
      al3N(3,4) = 2758881020544.d0/386414215625.d0
      al3N(4,4) = -469180684224.d0/386414215625.d0
      al3N(1,5) = 14032923046.d0/9273941175.d0
      al3N(2,5) = -23625504563402.d0/5332516175625.d0
      al3N(3,5) = 39862733245824.d0/8887526959375.d0
      al3N(4,5) = -605874357648.d0/386414215625.d0
      al3N(1,6) = 714229608.d0/3091313725.d0
      al3N(2,6) = -1625924907096.d0/1777505391875.d0
      al3N(3,6) = 11450529075456.d0/8887526959375.d0
      al3N(4,6) = -233665854912.d0/386414215625.d0
    
    
      al3D(1,1) = 1.d0/1.d0
      al3D(2,1) = 1.d0/1.d0
      al3D(3,1) = 1.d0/1.d0
      al3D(4,1) = 1.d0/1.d0
      al3D(5,1) = 1.d0/1.d0
      al3D(1,2) = 252019032.d0/123652549.d0
      al3D(2,2) = 252019032.d0/123652549.d0
      al3D(3,2) = 252019032.d0/123652549.d0
      al3D(4,2) = 252019032.d0/123652549.d0
      al3D(5,2) = 252019032.d0/123652549.d0
      al3D(1,3) = 151158380.d0/123652549.d0
      al3D(2,3) = 151158380.d0/123652549.d0
      al3D(3,3) = 151158380.d0/123652549.d0
      al3D(4,3) = 151158380.d0/123652549.d0
      al3D(5,3) = 151158380.d0/123652549.d0
    
    
c   This secondary predictor uses information from the new and old step 
    
    
c   U^{(n+1,4)} = al4_{1}*U^{(n  ,4)} + al4_{2}*U^{(n  ,5)} +
c               + al4_{3}*U^{(n  ,6)} + al4_{4}*U^{(n+1,2)} +
c               + al4_{5}*U^{(n+1,3)}
    
    
c   al4_{i} = \sum_{j=1}^{2*(order)} al4N(i,j)*r^{j} / 
c             \sum_{j=1}^{2*(order)} al4D(i,j)*r^{j} 
    
    
      al4N(3,1) = 17373884869672.d0/5111972803487.d0
      al4N(4,1) = 50855473860817.d0/10956348499794.d0
      al4N(5,1) = -117184182747644.d0/16644748100051.d0
      al4N(1,2) = -5438750468759.d0/3660786545912.d0
      al4N(2,2) = -3819405519603.d0/7956768730603.d0
      al4N(3,2) = 43297578418502.d0/14321563369351.d0
      al4N(4,2) = 22787025354911.d0/8976368084534.d0
      al4N(5,2) = -100903679058503.d0/46864611777020.d0
      al4N(1,3) = 961636147200.d0/1646934976121.d0
      al4N(2,3) = -6171567095808.d0/1646934976121.d0
      al4N(3,3) = 19971769123009.d0/7727660946277.d0
      al4N(4,3) = 8532165762273.d0/7499476144339.d0
      al4N(5,3) = 2396085749927.d0/4968281853644.d0
      al4N(1,4) = 14600964631732.d0/8810131861307.d0
      al4N(2,4) = -25141679652759.d0/6904059682726.d0
      al4N(3,4) = 18076449439006.d0/14809957156917.d0
      al4N(4,4) = 4106795411485.d0/7453642498521.d0
      al4N(5,4) = 5394407214677.d0/6737355498842.d0
      al4N(1,5) = 4286714975627.d0/10388638671097.d0
      al4N(2,5) = 1127715484801.d0/8458615186770.d0
      al4N(3,5) = -7223635457121.d0/9713562311645.d0
      al4N(4,5) = 1059591699648.d0/12531026992225.d0
      al4N(5,5) = 101712407400.d0/501241079689.d0
      al4N(1,6) = -152608335120.d0/1646934976121.d0
      al4N(2,6) = 4897031909184.d0/8234674880605.d0
      al4N(3,6) = -179738705808.d0/358029342635.d0
      al4N(1,7) = -169564816800.d0/1646934976121.d0
      al4N(2,7) = 429564202560.d0/1646934976121.d0
      al4N(3,7) = -11304321120.d0/71605868527.d0
    
    
      al4D(1,1) = 1.d0/1.d0
      al4D(2,1) = 1.d0/1.d0
      al4D(3,1) = 1.d0/1.d0
      al4D(4,1) = 1.d0/1.d0
      al4D(5,1) = 1.d0/1.d0
      al4D(1,2) = 103328303256.d0/71605868527.d0
      al4D(2,2) = 103328303256.d0/71605868527.d0
      al4D(3,2) = 103328303256.d0/71605868527.d0
      al4D(4,2) = 103328303256.d0/71605868527.d0
      al4D(5,2) = 103328303256.d0/71605868527.d0
      al4D(1,3) = 74543135229.d0/71605868527.d0
      al4D(2,3) = 74543135229.d0/71605868527.d0
      al4D(3,3) = 74543135229.d0/71605868527.d0
      al4D(4,3) = 74543135229.d0/71605868527.d0
      al4D(5,3) = 74543135229.d0/71605868527.d0
      al4D(1,4) = 42098769138.d0/71605868527.d0
      al4D(2,4) = 42098769138.d0/71605868527.d0
      al4D(3,4) = 42098769138.d0/71605868527.d0
      al4D(4,4) = 42098769138.d0/71605868527.d0
      al4D(5,4) = 42098769138.d0/71605868527.d0
      al4D(1,5) = 6428066472.d0/71605868527.d0
      al4D(2,5) = 6428066472.d0/71605868527.d0
      al4D(3,5) = 6428066472.d0/71605868527.d0
      al4D(4,5) = 6428066472.d0/71605868527.d0
      al4D(5,5) = 6428066472.d0/71605868527.d0
    
c   Dense output 
    
c   bd_{i} = \sum_{j=1}^{dense_order}  bdense_{ij}*x^{j} 
    
c   U[t^{n} + x*(Delta t)] = U^{(n)}  
c           + (Delta t)\sum_{i=1}^{s} bd_{i}*F^{(i)} 
    
    
      bD(1,1) = 2730611.d0/1312230.d0
      bD(1,2) = -303504119.d0/29218988.d0
      bD(1,3) = 1661475152.d0/109571205.d0
      bD(1,4) = -1468075291.d0/219142410.d0
      bD(2,1) = 1700.d0/189.d0
      bD(2,2) = -218347.d0/3507.d0
      bD(2,3) = 1026182.d0/10521.d0
      bD(2,4) = -1397323.d0/31563.d0
      bD(3,1) = -409375.d0/53784.d0
      bD(3,2) = 38609375.d0/665328.d0
      bD(3,3) = -242271875.d0/2619729.d0
      bD(3,4) = 2656240625.d0/62873496.d0
      bD(4,1) = 744925.d0/462024.d0
      bD(4,2) = -27907375.d0/1905136.d0
      bD(4,3) = 87388925.d0/3214917.d0
      bD(4,4) = -1038976925.d0/77158008.d0
      bD(5,1) = -60568.d0/8211.d0
      bD(5,2) = 23827180.d0/457079.d0
      bD(5,3) = -114128192.d0/1371237.d0
      bD(5,4) = 52384088.d0/1371237.d0
      bD(6,1) = 33.d0/10.d0
      bD(6,2) = -15273.d0/668.d0
      bD(6,3) = 29916.d0/835.d0
      bD(6,4) = -26743.d0/1670.d0
   

      elseif(icase.eq.2) then     !    The scheme you're using is: QESDIRK436L2SA
    
      ns = 6 

c     The scheme coefficients, A, b, bh, and c  are given by
    
      ai(2,1) = 8.d0/75.d0
      ai(2,2) = 8.d0/75.d0
      ai(3,1) = 605497861978.d0/9136257189845.d0
      ai(3,2) = -2127848798551.d0/10702252975294.d0
      ai(3,3) = 8.d0/25.d0
      ai(4,1) = -3005106686955.d0/6150333508049.d0
      ai(4,2) = -68662668193799.d0/11091168490809.d0
      ai(4,3) = 80786898110822.d0/11737001380747.d0
      ai(4,4) = 8.d0/25.d0
      ai(5,1) = -26162805558846.d0/8363194173203.d0
      ai(5,2) = -291987295964487.d0/9066074244437.d0
      ai(5,3) = 384682892278670.d0/10959450712301.d0
      ai(5,4) = 13555548451102.d0/14148104892819.d0
      ai(5,5) = 8.d0/25.d0
      ai(6,1) = 540088238697.d0/4693226184761.d0
      ai(6,3) = 1094762490994.d0/7880993776667.d0
      ai(6,4) = 4016564763003.d0/7185357966874.d0
      ai(6,5) = -411820258827.d0/3096789411938.d0
      ai(6,6) = 8.d0/25.d0
    
      be(1) = 540088238697.d0/4693226184761.d0
      be(2) = 0.d0/1.d0
      be(3) = 1094762490994.d0/7880993776667.d0
      be(4) = 4016564763003.d0/7185357966874.d0
      be(5) = -411820258827.d0/3096789411938.d0
      be(6) = 8.d0/25.d0
    
      beh(1) = -374484326677.d0/8442488809460.d0
      beh(2) = -41125091159938.d0/25963879779069.d0
      beh(3) = 24025088270494.d0/12927594097169.d0
      beh(4) = 5193917461301.d0/8985383982321.d0
      beh(5) = -1843122001830.d0/16078617943063.d0
      beh(6) = 2439572212972.d0/7960792257433.d0
    
      bi(1) = 540088238697.d0/4693226184761.d0
      bi(2) = 0.d0/1.d0
      bi(3) = 1094762490994.d0/7880993776667.d0
      bi(4) = 4016564763003.d0/7185357966874.d0
      bi(5) = -411820258827.d0/3096789411938.d0
      bi(6) = 8.d0/25.d0
    
      bih(1) = -374484326677.d0/8442488809460.d0
      bih(2) = -41125091159938.d0/25963879779069.d0
      bih(3) = 24025088270494.d0/12927594097169.d0
      bih(4) = 5193917461301.d0/8985383982321.d0
      bih(5) = -1843122001830.d0/16078617943063.d0
      bih(6) = 2439572212972.d0/7960792257433.d0
    
c   The primary predictor is now given in two parts; eta and B
    
c   U^{(n+1,i+1)} = \sum_{j=0}^{order}eta_{ij}*r^{j} +  
c                 + \sum_{j=1}^{S-1}\sum_{k=0}^{order}B_{ijk}*r^{k}U^{(n,j+1)} 
    
      svpB(2,6,0) = 1.d0/1.d0
      svpB(3,6,0) = 1.d0/1.d0
      svpB(4,6,0) = 1.d0/1.d0
      svpB(5,6,0) = 1.d0/1.d0
      svpB(6,6,0) = 1.d0/1.d0
      svpB(2,2,1) = -7434808665603555.d0/9738477798251.d0
      svpB(3,2,1) = -5539058360398754.d0/11618548485105.d0
      svpB(4,2,1) = 39654491347958446.d0/11555768347435.d0
      svpB(5,2,1) = 181462141912838692.d0/8221526358953.d0
      svpB(6,2,1) = -11957964214393759.d0/13948859216214.d0
      svpB(2,3,1) = 7611586547869294.d0/8920594527397.d0
      svpB(3,3,1) = 4386479054425647.d0/8225300248102.d0
      svpB(4,3,1) = -60273901707948946.d0/15767762355131.d0
      svpB(5,3,1) = -250828456017122087.d0/10193580961079.d0
      svpB(6,3,1) = 15487796373272975.d0/16055836356176.d0
      svpB(2,4,1) = -57064021328902.d0/6196250797211.d0
      svpB(3,4,1) = -29993787314137.d0/4810389049043.d0
      svpB(4,4,1) = 115584742963983.d0/4105160004940.d0
      svpB(5,4,1) = 1781550979297471.d0/8825842411462.d0
      svpB(6,4,1) = -171750077361643.d0/10025328486495.d0
      svpB(2,5,1) = -86311485209687.d0/16152641864906.d0
      svpB(3,5,1) = -21317875757559.d0/6561934038245.d0
      svpB(4,5,1) = 61880507411516.d0/2341635374587.d0
      svpB(5,5,1) = 2771428086915300.d0/16675759871803.d0
      svpB(6,5,1) = -35806796669701.d0/7514379047075.d0
      svpB(2,6,1) = 115415483495444.d0/8525531159265.d0
      svpB(3,6,1) = 119711886406591.d0/13947918979868.d0
      svpB(4,6,1) = -1230765466413885.d0/21471294570181.d0
      svpB(5,6,1) = -2373727781993353.d0/6343157133377.d0
      svpB(6,6,1) = 148743169983109.d0/8743966818208.d0
      svpB(2,2,2) = 1194399771414241.d0/14700698573962.d0
      svpB(3,2,2) = 753698350794587.d0/12015025118895.d0
      svpB(4,2,2) = 12748661878462255.d0/26034375650777.d0
      svpB(5,2,2) = 18536096487102889.d0/9530501990468.d0
      svpB(6,2,2) = 26196850674773025.d0/14674211838863.d0
      svpB(2,3,2) = -767702347954113.d0/8586424496552.d0
      svpB(3,3,2) = -354620662139171.d0/5137151571771.d0
      svpB(4,3,2) = -2734078147582287.d0/5073699167078.d0
      svpB(5,3,2) = -27400740183287666.d0/12802388877593.d0
      svpB(6,3,2) = -15150838485197579.d0/7712116713794.d0
      svpB(2,4,2) = -515216514758.d0/5318352224079.d0
      svpB(3,4,2) = -361357979957.d0/4831298017400.d0
      svpB(4,4,2) = -5514350294843.d0/9444434842985.d0
      svpB(5,4,2) = -5549954513793.d0/2393238487519.d0
      svpB(6,4,2) = -14868661029533.d0/6985164210555.d0
      svpB(2,5,2) = 5377865767535.d0/6404795703259.d0
      svpB(3,5,2) = 12190849623797.d0/18804775568152.d0
      svpB(4,5,2) = 31800490558079.d0/6283815799779.d0
      svpB(5,5,2) = 403338453425501.d0/20066618120420.d0
      svpB(6,5,2) = 156038578502294.d0/8457553374701.d0
      svpB(2,6,2) = -19758784721423.d0/14126458386748.d0
      svpB(3,6,2) = -18034730297452.d0/16700196218351.d0
      svpB(4,6,2) = -38101463517843.d0/4519693929055.d0
      svpB(5,6,2) = -83036143816563.d0/2479987381910.d0
      svpB(6,6,2) = -630665662207213.d0/20520582718778.d0
      svpB(2,2,3) = 545817172836139.d0/9216734577429.d0
      svpB(3,2,3) = 59995039945747.d0/2239983279613.d0
      svpB(4,2,3) = 7312250842258879.d0/12517387012214.d0
      svpB(5,2,3) = 66624196939920562.d0/14408440659441.d0
      svpB(6,2,3) = 156656345444247079.d0/38525234689675.d0
      svpB(2,3,3) = -247042989971029.d0/3754587518749.d0
      svpB(3,3,3) = -334401525240253.d0/11237183974790.d0
      svpB(4,3,3) = -37823595653459917.d0/58275394570788.d0
      svpB(5,3,3) = -137119304875123951.d0/26689682638995.d0
      svpB(6,3,3) = -42407458944746043.d0/9386406627227.d0
      svpB(2,4,3) = 5797795397253.d0/13497520191809.d0
      svpB(3,4,3) = 1987307028959.d0/10229528430650.d0
      svpB(4,4,3) = 38936250872486.d0/9189205114011.d0
      svpB(5,4,3) = 69741516307817.d0/2079397875082.d0
      svpB(6,4,3) = 497148837255801.d0/16855628295982.d0
      svpB(2,5,3) = 3595109637688.d0/6916875640895.d0
      svpB(3,5,3) = 3096953975294.d0/13174438784747.d0
      svpB(4,5,3) = 14149547275467.d0/2759771250859.d0
      svpB(5,5,3) = 679029118144177.d0/16731748809367.d0
      svpB(6,5,3) = 135112530985133.d0/3785827961319.d0
      svpB(2,6,3) = -5076327536726.d0/4756458745451.d0
      svpB(3,6,3) = -10351531005983.d0/21445601414240.d0
      svpB(4,6,3) = -298841779125428.d0/28386256644975.d0
      svpB(5,6,3) = -1625423025710692.d0/19505433895813.d0
      svpB(6,6,3) = -794210673234819.d0/10837702725028.d0
      svpB(2,2,4) = 507383452403427.d0/39242209927450.d0
      svpB(3,2,4) = -49358508221957.d0/9425956918773.d0
      svpB(4,2,4) = 588247782390826.d0/8809817872959.d0
      svpB(5,2,4) = 28181786054500141.d0/10619570363121.d0
      svpB(6,2,4) = 15122285170888414.d0/6758780975377.d0
      svpB(2,3,4) = -48649057486955.d0/3372950623257.d0
      svpB(3,3,4) = 51694727079137.d0/8852815548266.d0
      svpB(4,3,4) = -1324541945635984.d0/17761112960819.d0
      svpB(5,3,4) = -82786541562609367.d0/27958341771278.d0
      svpB(6,3,4) = -9605370771923313.d0/3848051885060.d0
      svpB(2,4,4) = 1506774188039.d0/11835900594831.d0
      svpB(3,4,4) = -385757060333.d0/7795438620081.d0
      svpB(4,4,4) = 3628099262584.d0/4855142983375.d0
      svpB(5,4,4) = 201347801941209.d0/7498984193443.d0
      svpB(6,4,4) = 117128255862632.d0/5258076818485.d0
      svpB(2,5,4) = 1213752668959.d0/11661792907042.d0
      svpB(3,5,4) = -157941011597.d0/4188763662512.d0
      svpB(4,5,4) = 7937600247683.d0/10871861101499.d0
      svpB(5,5,4) = 195491766616426.d0/8533881598023.d0
      svpB(6,5,4) = 259868768387891.d0/14018105378600.d0
      svpB(2,6,4) = -1472959851030.d0/6405868195373.d0
      svpB(3,6,4) = 421936802368.d0/4806485976303.d0
      svpB(4,6,4) = -17675638540147.d0/12458067239952.d0
      svpB(5,6,4) = -579038683193116.d0/11804839145765.d0
      svpB(6,6,4) = -314421512020375.d0/7778113839049.d0
    
c   U^{(n+1,i+1)} = \sum_{j=0}^{order}eta_{ij}*r^{j} +  
c                 + \sum_{j=1}^{S-1}\sum_{k=0}^{order}B_{ijk}*r^{k}U^{(n,j+1)} 
    
      svpB(2,1,1) = -2253762174705014.d0/25380773384213.d0
      svpB(3,1,1) = -144269233688933.d0/2592570023600.d0
      svpB(4,1,1) = 4015813215651479.d0/10198445767636.d0
      svpB(5,1,1) = 25001581082506343.d0/9838920180710.d0
      svpB(6,1,1) = -547342086122669.d0/5341832612639.d0
      svpB(2,1,2) = 7510929898799.d0/851874192014.d0
      svpB(3,1,2) = 42721587169882.d0/6275784013127.d0
      svpB(4,1,2) = 530854733407941.d0/9989683979759.d0
      svpB(5,1,2) = 5279340568304789.d0/25013281197391.d0
      svpB(6,1,2) = 724146283486493.d0/3737881226306.d0
      svpB(2,1,3) = 80084162641606.d0/11961162511259.d0
      svpB(3,1,3) = 32855125524652.d0/10849992231023.d0
      svpB(4,1,3) = 550393700635195.d0/8333605161713.d0
      svpB(5,1,3) = 4615361101895577.d0/8828512722595.d0
      svpB(6,1,3) = 1418543407728454.d0/3085579913665.d0
      svpB(2,1,4) = 11676342553607.d0/7824315420637.d0
      svpB(3,1,4) = -7545121298075.d0/12502209086062.d0
      svpB(4,1,4) = 64969581690837.d0/8388601447769.d0
      svpB(5,1,4) = 3517879424032543.d0/11473789495542.d0
      svpB(6,1,4) = 2657521186584219.d0/10286651190638.d0
    
c  The secondary predictor is used well into the next step
c  It uses only information from the new step 
    
c   Z_{0}^{(n+1,i)} = \sum_{j=2}^{i-1} \alpha_{ij} Z_{0}^{(n+1,j)} 
    
      alpha(3,2) = 6422274662979.d0/7309005751876.d0
      alpha(4,2) = 120276492383629.d0/3770638964174.d0
      alpha(4,3) = -1083533752239991.d0/32336165880218.d0
      alpha(5,2) = 1097129296965975.d0/7413256282517.d0
      alpha(5,3) = -706969143941492.d0/4308757153593.d0
      alpha(5,4) = 7890335017315.d0/18130153042293.d0
      alpha(6,2) = 38438471859897.d0/839018322652.d0
      alpha(6,3) = -401969876064773.d0/7854340830358.d0
      alpha(6,4) = 3702428761097.d0/7078989906133.d0
      alpha(6,5) = 1396203549671.d0/2669484771837.d0
    
c   This secondary predictor uses information from the new and old step 
    
    
c   U^{(n+1,3)} = al3_{1}*U^{(n  ,3)} + al3_{2}*U^{(n  ,4)} +
c               + al3_{3}*U^{(n  ,5)} + al3_{4}*U^{(n  ,6)} +
c               + al3_{5}*U^{(n+1,2)}
    
    
c   al3_{i} = \sum_{j=1}^{2*(order)} al3N(i,j)*r^{j} / 
c             \sum_{j=1}^{2*(order)} al3D(i,j)*r^{j} 
    
    
      al3N(1,1) = -44581474551.d0/9445280557853.d0
      al3N(2,1) = 125791452133.d0/5560416844349.d0
      al3N(3,1) = -12893224131922.d0/12256564736041.d0
      al3N(4,1) = 6112950699321.d0/4843998733373.d0
      al3N(5,1) = 4648812485185.d0/6021169952558.d0
      al3N(1,2) = 362416114790.d0/14258407860089.d0
      al3N(2,2) = -1101419494510.d0/8705367776147.d0
      al3N(3,2) = -6658861458763.d0/17283095181369.d0
      al3N(4,2) = 6470897335788.d0/10014300353437.d0
      al3N(5,2) = 2125686186628.d0/16111695750689.d0
      al3N(1,3) = 229048686371.d0/9081797724800.d0
      al3N(2,3) = -526363799557.d0/6967768342385.d0
      al3N(3,3) = 16859141458357.d0/66772161228200.d0
      al3N(4,3) = -1825169416991.d0/10084661233721.d0
      al3N(5,3) = -35480733715.d0/4978213646872.d0
      al3N(1,4) = 454223079392.d0/36550379694679.d0
      al3N(2,4) = -29831970162.d0/1612191109637.d0
      al3N(3,4) = 541386629890.d0/4355933011457.d0
      al3N(4,4) = -434692760047.d0/3677276948858.d0
      al3N(1,5) = -22820657357.d0/15089762993135.d0
      al3N(2,5) = 77937044077.d0/10319893792534.d0
      al3N(3,5) = 384708024509.d0/7110737589269.d0
      al3N(4,5) = -402030428816.d0/6684664765077.d0
      al3N(1,6) = -12660715562.d0/15086347745225.d0
      al3N(2,6) = 2770568593.d0/1175143931830.d0
      al3N(3,6) = 69431568569.d0/6891585945120.d0
      al3N(4,6) = -110273281329.d0/9511846623400.d0
    
    
      al3D(1,1) = 1.d0/1.d0
      al3D(2,1) = 1.d0/1.d0
      al3D(3,1) = 1.d0/1.d0
      al3D(4,1) = 1.d0/1.d0
      al3D(5,1) = 1.d0/1.d0
      al3D(1,2) = 4099692434849.d0/14053802107402.d0
      al3D(2,2) = 4099692434849.d0/14053802107402.d0
      al3D(3,2) = 4099692434849.d0/14053802107402.d0
      al3D(4,2) = 4099692434849.d0/14053802107402.d0
      al3D(5,2) = 4099692434849.d0/14053802107402.d0
      al3D(1,3) = 187246148419.d0/13323673727852.d0
      al3D(2,3) = 187246148419.d0/13323673727852.d0
      al3D(3,3) = 187246148419.d0/13323673727852.d0
      al3D(4,3) = 187246148419.d0/13323673727852.d0
      al3D(5,3) = 187246148419.d0/13323673727852.d0
    
    
c   This secondary predictor uses information from the new and old step 
    
    
c   U^{(n+1,4)} = al4_{1}*U^{(n  ,4)} + al4_{2}*U^{(n  ,5)} +
c               + al4_{3}*U^{(n  ,6)} + al4_{4}*U^{(n+1,2)} +
c               + al4_{5}*U^{(n+1,3)}
    
    
c   al4_{i} = \sum_{j=1}^{2*(order)} al4N(i,j)*r^{j} / 
c             \sum_{j=1}^{2*(order)} al4D(i,j)*r^{j} 
    
    
      al4N(1,1) = -1.d0/653155474825006662015530.d0
      al4N(2,1) = 1.d0/5516180969363842963578.d0
      al4N(3,1) = -114151788274401.d0/12988716289984.d0
      al4N(4,1) = -51227372312417.d0/7629032402710.d0
      al4N(5,1) = 135311488700123.d0/8199042610759.d0
      al4N(1,2) = -149029087154.d0/7326886952115.d0
      al4N(2,2) = -368579054827348.d0/14064101845963.d0
      al4N(3,2) = 194697263628951.d0/16250031663619.d0
      al4N(4,2) = 14740196254463.d0/2787885197273.d0
      al4N(5,2) = 107093784389302.d0/10762017178717.d0
      al4N(1,3) = 6766162452308.d0/15478547960371.d0
      al4N(2,3) = -2151039686647324.d0/41558336351801.d0
      al4N(3,3) = 1394631515805871.d0/27775060497112.d0
      al4N(4,3) = 156727580121139.d0/21548121295581.d0
      al4N(5,3) = -99102727489949.d0/17467016334508.d0
      al4N(1,4) = -6447509522833.d0/5185296146208.d0
      al4N(2,4) = -225236847004048.d0/11960396507019.d0
      al4N(3,4) = 53897106457434.d0/2638074778397.d0
      al4N(4,4) = -771985497700.d0/4575557405269.d0
      al4N(5,4) = -1359050223275.d0/5525954701979.d0
      al4N(1,5) = 4928541700.d0/5410753950327.d0
      al4N(2,5) = 13034796537767.d0/11106481755551.d0
      al4N(3,5) = -3610746973390.d0/5086565168191.d0
      al4N(4,5) = -786962834248.d0/4284734667429.d0
      al4N(5,5) = -452556243263.d0/1441159405189.d0
      al4N(1,6) = -151330320449.d0/10604294275499.d0
      al4N(2,6) = 22023036562042.d0/13033316848511.d0
      al4N(3,6) = -12284800295940.d0/7332114449519.d0
      al4N(1,7) = 683969414842.d0/16388612949169.d0
      al4N(2,7) = 4572008404247.d0/10067543248919.d0
      al4N(3,7) = -2567708973447.d0/5178211671202.d0
    
    
      al4D(1,1) = 1.d0/1.d0
      al4D(2,1) = 1.d0/1.d0
      al4D(3,1) = 1.d0/1.d0
      al4D(4,1) = 1.d0/1.d0
      al4D(5,1) = 1.d0/1.d0
      al4D(1,2) = 11697759734512.d0/11789168721383.d0
      al4D(2,2) = 11697759734512.d0/11789168721383.d0
      al4D(3,2) = 11697759734512.d0/11789168721383.d0
      al4D(4,2) = 11697759734512.d0/11789168721383.d0
      al4D(5,2) = 11697759734512.d0/11789168721383.d0
      al4D(1,3) = 14800162518005.d0/30270921842057.d0
      al4D(2,3) = 14800162518005.d0/30270921842057.d0
      al4D(3,3) = 14800162518005.d0/30270921842057.d0
      al4D(4,3) = 14800162518005.d0/30270921842057.d0
      al4D(5,3) = 14800162518005.d0/30270921842057.d0
      al4D(1,4) = -362292843013.d0/6089063314218.d0
      al4D(2,4) = -362292843013.d0/6089063314218.d0
      al4D(3,4) = -362292843013.d0/6089063314218.d0
      al4D(4,4) = -362292843013.d0/6089063314218.d0
      al4D(5,4) = -362292843013.d0/6089063314218.d0
      al4D(1,5) = -806525308961.d0/24427573525771.d0
      al4D(2,5) = -806525308961.d0/24427573525771.d0
      al4D(3,5) = -806525308961.d0/24427573525771.d0
      al4D(4,5) = -806525308961.d0/24427573525771.d0
      al4D(5,5) = -806525308961.d0/24427573525771.d0
    
c   Dense output 
    
c   bd_{i} = \sum_{j=1}^{dense_order}  bdense_{ij}*x^{j} 
    
c   U[t^{n} + x*(Delta t)] = U^{(n)}  
c           + (Delta t)\sum_{i=1}^{s} bd_{i}*F^{(i)} 
    
    
      bD(1,1) = 9045775731155.d0/13719517693781.d0
      bD(1,2) = -5372340945483.d0/4213705305623.d0
      bD(1,3) = 14404239339754.d0/13954339267035.d0
      bD(1,4) = -846216657906.d0/2806404445669.d0
      bD(2,1) = -1.d0/7052714834817224339197.d0
      bD(2,2) = 1.d0/854113406724811764394.d0
      bD(2,3) = -1.d0/520257875360813648575.d0
      bD(2,4) = 1.d0/1119979724547366783262.d0
      bD(3,1) = 5139185493805.d0/6457161246403.d0
      bD(3,2) = -9665876204501.d0/6280527642738.d0
      bD(3,3) = 1365936830285.d0/1096235530947.d0
      bD(3,4) = -4021622342545.d0/11049033357932.d0
      bD(4,1) = -8303389601397.d0/10375793610670.d0
      bD(4,2) = 69395979835253.d0/12457122083915.d0
      bD(4,3) = -62796494422246.d0/9653860459083.d0
      bD(4,4) = 48120001405582.d0/20983069448621.d0
      bD(5,1) = -20549818646764.d0/11975853676037.d0
      bD(5,2) = 147398745884552.d0/10278279541867.d0
      bD(5,3) = -274948833053050.d0/11424917200699.d0
      bD(5,4) = 715533616246667.d0/63277471806833.d0
      bD(6,1) = 8517941425826.d0/4132961650537.d0
      bD(6,2) = -138534237134522.d0/8102556638089.d0
      bD(6,3) = 156138389025553.d0/5518767234093.d0
      bD(6,4) = -60416963794213.d0/4670581570396.d0
    

      elseif(icase.eq.3) then     !    The scheme you're using is: ESDIRK436L2SA
    
      ns = 6

c   The scheme coefficients, A, b, bh, and c  are given by
    
      ai(2,1) = 1.d0/4.d0
      ai(2,2) = 1.d0/4.d0
      ai(3,1) = -1356991263433.d0/26208533697614.d0
      ai(3,2) = -1356991263433.d0/26208533697614.d0
      ai(3,3) = 1.d0/4.d0
      ai(4,1) = -162708388469.d0/2125389860943.d0
      ai(4,2) = -162708388469.d0/2125389860943.d0
      ai(4,3) = 1768468916152.d0/3348680272939.d0
      ai(4,4) = 1.d0/4.d0
      ai(5,1) = -1410042975763.d0/1938452943086.d0
      ai(5,2) = -1410042975763.d0/1938452943086.d0
      ai(5,3) = 20915361222208.d0/13195852609937.d0
      ai(5,4) = 7240912463611.d0/10974111771892.d0
      ai(5,5) = 1.d0/4.d0
      ai(6,1) = -140404485182.d0/9007427031765.d0
      ai(6,2) = -140404485182.d0/9007427031765.d0
      ai(6,3) = 1402990318619.d0/3619147572429.d0
      ai(6,4) = 1394850835967.d0/2779846451479.d0
      ai(6,5) = -990969975725.d0/9154032505244.d0
      ai(6,6) = 1.d0/4.d0
    
      be(1) = -140404485182.d0/9007427031765.d0
      be(2) = -140404485182.d0/9007427031765.d0
      be(3) = 1402990318619.d0/3619147572429.d0
      be(4) = 1394850835967.d0/2779846451479.d0
      be(5) = -990969975725.d0/9154032505244.d0
      be(6) = 1.d0/4.d0
    
      beh(1) = -480923228411.d0/4982971448372.d0
      beh(2) = -480923228411.d0/4982971448372.d0
      beh(3) = 6709447293961.d0/12833189095359.d0
      beh(4) = 3513175791894.d0/6748737351361.d0
      beh(5) = -498863281070.d0/6042575550617.d0
      beh(6) = 2077005547802.d0/8945017530137.d0
    
      bi(1) = -140404485182.d0/9007427031765.d0
      bi(2) = -140404485182.d0/9007427031765.d0
      bi(3) = 1402990318619.d0/3619147572429.d0
      bi(4) = 1394850835967.d0/2779846451479.d0
      bi(5) = -990969975725.d0/9154032505244.d0
      bi(6) = 1.d0/4.d0
    
      bih(1) = -480923228411.d0/4982971448372.d0
      bih(2) = -480923228411.d0/4982971448372.d0
      bih(3) = 6709447293961.d0/12833189095359.d0
      bih(4) = 3513175791894.d0/6748737351361.d0
      bih(5) = -498863281070.d0/6042575550617.d0
      bih(6) = 2077005547802.d0/8945017530137.d0
    
c   The primary predictor is now given in two parts; eta and B
    
c   U^{(n+1,i+1)} = \sum_{j=0}^{order}eta_{ij}*r^{j} +  
c                 + \sum_{j=1}^{S-1}\sum_{k=0}^{order}B_{ijk}*r^{k}U^{(n,j+1)} 
    
      svpB(2,6,0) = 1.d0/1.d0
      svpB(3,6,0) = 1.d0/1.d0
      svpB(4,6,0) = 1.d0/1.d0
      svpB(5,6,0) = 1.d0/1.d0
      svpB(6,6,0) = 1.d0/1.d0
      svpB(2,2,1) = -543586075368876.d0/24533375850871.d0
      svpB(3,2,1) = 62624216366346.d0/11522505303659.d0
      svpB(4,2,1) = 73114831274098.d0/7763092679493.d0
      svpB(5,2,1) = 218176496509409.d0/2992332624706.d0
      svpB(6,2,1) = 34438766894071.d0/7069672812925.d0
      svpB(2,3,1) = 1354170239561327.d0/13125859662734.d0
      svpB(3,3,1) = -222155604493141.d0/11687035997772.d0
      svpB(4,3,1) = -196221583660269.d0/8090759425333.d0
      svpB(5,3,1) = -3077401735686405.d0/11123867448502.d0
      svpB(6,3,1) = 35105511941288.d0/10657885425011.d0
      svpB(2,4,1) = -325160135838757.d0/9203294859050.d0
      svpB(3,4,1) = 14188317307243.d0/3404671390812.d0
      svpB(4,4,1) = 12170310062455.d0/11991146759643.d0
      svpB(5,4,1) = 1959300999917873.d0/27454498151040.d0
      svpB(6,4,1) = -476349274131952.d0/44148537643513.d0
      svpB(2,5,1) = -84563175949208.d0/8298692591421.d0
      svpB(3,5,1) = 7652624134515.d0/3008786236774.d0
      svpB(4,5,1) = 25531672694023.d0/5714218858866.d0
      svpB(5,5,1) = 79982206901732.d0/2354486829463.d0
      svpB(6,5,1) = 29831696078735.d0/12319884342642.d0
      svpB(2,6,1) = 213365900458718.d0/7319787752343.d0
      svpB(3,6,1) = -50562339094378.d0/10038194450979.d0
      svpB(4,6,1) = -59896448544569.d0/10302845468187.d0
      svpB(5,6,1) = -1758009271570079.d0/23492126061027.d0
      svpB(6,6,1) = 15910435145759.d0/6895893816188.d0
      svpB(2,2,2) = 95774543217135.d0/11845332587134.d0
      svpB(3,2,2) = 3947717417399.d0/5691471872288.d0
      svpB(4,2,2) = 226459886469986.d0/17925382392605.d0
      svpB(5,2,2) = 212358350081899.d0/6070714966671.d0
      svpB(6,2,2) = 352984813179585.d0/10914232451206.d0
      svpB(2,3,2) = -487370910666614.d0/29820216034581.d0
      svpB(3,3,2) = -12805294925308.d0/9133182438149.d0
      svpB(4,3,2) = -34839761891825.d0/1364288861359.d0
      svpB(5,3,2) = -710730346839796.d0/10051465136591.d0
      svpB(6,3,2) = -1333972089718493.d0/20405062668337.d0
      svpB(2,4,2) = -1672027644349.d0/8418424476640.d0
      svpB(3,4,2) = -100182514186.d0/5879774351901.d0
      svpB(4,4,2) = -13931442659854.d0/44891429216831.d0
      svpB(5,4,2) = -29832048246271.d0/34717126297731.d0
      svpB(6,4,2) = -5972621778671.d0/7517828060620.d0
      svpB(2,5,2) = 36542625360245.d0/8359464793717.d0
      svpB(3,5,2) = 2438834123657.d0/6503428304942.d0
      svpB(4,5,2) = 42923316155063.d0/6284228510731.d0
      svpB(5,5,2) = 212747675079367.d0/11249076764099.d0
      svpB(6,5,2) = 119772524689303.d0/6849769231989.d0
      svpB(2,6,2) = -34345469493625.d0/5656954970946.d0
      svpB(3,6,2) = -5716313309255.d0/10975154812546.d0
      svpB(4,6,2) = -132709886513449.d0/13989317149995.d0
      svpB(5,6,2) = -271631680370333.d0/10341106592324.d0
      svpB(6,6,2) = -69735957481120.d0/2871508070371.d0
      svpB(2,2,3) = 105460779055436.d0/10040533249767.d0
      svpB(3,2,3) = -51230208233919.d0/40202984374502.d0
      svpB(4,2,3) = 360332848995098.d0/23865020160581.d0
      svpB(5,2,3) = 2581042037877120.d0/40386351410773.d0
      svpB(6,2,3) = 452603906003813.d0/8079516517730.d0
      svpB(2,3,3) = -173347480432044.d0/4740940012267.d0
      svpB(3,3,3) = 38872106800274.d0/8762970996679.d0
      svpB(4,3,3) = -630590434239034.d0/11997381017437.d0
      svpB(5,3,3) = -1993290904636691.d0/8959656297264.d0
      svpB(6,3,3) = -516344692216540.d0/2647816691671.d0
      svpB(2,4,3) = 177634027729091.d0/18914581702741.d0
      svpB(3,4,3) = -12669099968651.d0/11119434722230.d0
      svpB(4,4,3) = 73412517732841.d0/5437920366413.d0
      svpB(5,4,3) = 235841377079295.d0/4127283663166.d0
      svpB(6,4,3) = 1798455292918204.d0/35906369557443.d0
      svpB(2,5,3) = 75948709350054.d0/12318218606453.d0
      svpB(3,5,3) = -15855696973733.d0/21197235010084.d0
      svpB(4,5,3) = 55211420466898.d0/6229433821203.d0
      svpB(5,5,3) = 399657667632513.d0/10653428536939.d0
      svpB(6,5,3) = 208688369980480.d0/6346393961241.d0
      svpB(2,6,3) = -65593310504859.d0/5385824404636.d0
      svpB(3,6,3) = 3863557721727.d0/2614847786026.d0
      svpB(4,6,3) = -237786638881004.d0/13582266085681.d0
      svpB(5,6,3) = -966095892924733.d0/13037268616886.d0
      svpB(6,6,3) = -793080512302261.d0/12209874219040.d0
      svpB(2,2,4) = 25140295195094.d0/5158848483379.d0
      svpB(3,2,4) = -10002431837453.d0/6487913551906.d0
      svpB(4,2,4) = 50184302703491.d0/10920718547202.d0
      svpB(5,2,4) = 1453223957733016.d0/44988813456467.d0
      svpB(6,2,4) = 385119319021954.d0/14231358406581.d0
      svpB(2,3,4) = -88086114936308.d0/4159254612279.d0
      svpB(3,3,4) = 17841269087227.d0/2720779331679.d0
      svpB(4,3,4) = -154608916304891.d0/7442475778334.d0
      svpB(5,3,4) = -1790584990640997.d0/12311014438393.d0
      svpB(6,3,4) = -1739929164243847.d0/14475191534963.d0
      svpB(2,4,4) = 79044442885627.d0/11976996968624.d0
      svpB(3,4,4) = -15156173643869.d0/7585035644420.d0
      svpB(4,4,4) = 75600012929957.d0/11235410804417.d0
      svpB(5,4,4) = 197186570333883.d0/4201407336625.d0
      svpB(6,4,4) = 217407925811671.d0/5679111504114.d0
      svpB(2,5,4) = 30828636314353.d0/10885133337008.d0
      svpB(3,5,4) = -7268408731663.d0/8874903943965.d0
      svpB(4,5,4) = 30225320153715.d0/9736067843054.d0
      svpB(5,5,4) = 288574864668635.d0/13416364528896.d0
      svpB(6,5,4) = 114428198107621.d0/6680208554623.d0
      svpB(2,6,4) = -17879356334807.d0/2791297421965.d0
      svpB(3,6,4) = 58277302863925.d0/30493488257734.d0
      svpB(4,6,4) = -86374340455505.d0/12911915586657.d0
      svpB(5,6,4) = -663521469261769.d0/14252771718234.d0
      svpB(6,6,4) = -76928088586673.d0/2042236736885.d0
    
c   U^{(n+1,i+1)} = \sum_{j=0}^{order}eta_{ij}*r^{j} +  
c                 + \sum_{j=1}^{S-1}\sum_{k=0}^{order}B_{ijk}*r^{k}U^{(n,j+1)} 
    
      svpB(2,1,1) = -958900105692004.d0/14834577708743.d0
      svpB(3,1,1) = 125978163515225.d0/10586370289279.d0
      svpB(4,1,1) = 102152378871531.d0/6736132748278.d0
      svpB(5,1,1) = 2131824007193704.d0/12305954348637.d0
      svpB(6,1,1) = -20546653886902.d0/9764796125937.d0
      svpB(2,1,2) = 97489231232661.d0/9598424558198.d0
      svpB(3,1,2) = 8894156638657.d0/10207736485204.d0
      svpB(4,1,2) = 121132111045329.d0/7632776927620.d0
      svpB(5,1,2) = 263432657021125.d0/5994959484798.d0
      svpB(6,1,2) = 283615807289128.d0/6980937523449.d0
      svpB(2,1,3) = 245009687964347.d0/10801760848934.d0
      svpB(3,1,3) = -29359378089575.d0/10669019794537.d0
      svpB(4,1,3) = 651646731640954.d0/19985529939923.d0
      svpB(5,1,3) = 1571199710714081.d0/11384564694934.d0
      svpB(6,1,3) = 1372500928694174.d0/11345541261963.d0
      svpB(2,1,4) = 131261491130918.d0/9885167816321.d0
      svpB(3,1,4) = -40271937951939.d0/9799250074495.d0
      svpB(4,1,4) = 85311020189483.d0/6544847978412.d0
      svpB(5,1,4) = 1603878366668195.d0/17575766161473.d0
      svpB(6,1,4) = 224599442930539.d0/2978909507021.d0
    
c  The secondary predictor is used well into the next step
c  It uses only information from the new step 
    
c   Z_{0}^{(n+1,i)} = \sum_{j=2}^{i-1} \alpha_{ij} Z_{0}^{(n+1,j)} 
    
      alpha(3,2) = 1513744654945.d0/5168247530883.d0
      alpha(4,2) = 18037384292422.d0/10660759696825.d0
      alpha(4,3) = -19707882338879.d0/13061235440667.d0
      alpha(5,2) = 36510702415949.d0/10838657953295.d0
      alpha(5,3) = -148089091958132.d0/16162417422911.d0
      alpha(5,4) = 8372411349895.d0/7501686571312.d0
      alpha(6,2) = 4652010428369.d0/15234697295817.d0
      alpha(6,3) = -11355430753181.d0/5492328479418.d0
      alpha(6,4) = 32949210879299.d0/37466343683880.d0
      alpha(6,5) = 5591347549239.d0/9684332989580.d0
    
c   This secondary predictor uses information from the new and old step 
    
    
c   U^{(n+1,3)} = al3_{1}*U^{(n  ,3)} + al3_{2}*U^{(n  ,4)} +
c               + al3_{3}*U^{(n  ,5)} + al3_{4}*U^{(n  ,6)} +
c               + al3_{5}*U^{(n+1,2)}
    
    
c   al3_{i} = \sum_{j=1}^{2*(order)} al3N(i,j)*r^{j} / 
c             \sum_{j=1}^{2*(order)} al3D(i,j)*r^{j} 
    
    
      al3N(1,1) = -2034754227286.d0/12123849673973.d0
      al3N(2,1) = 10298564741129.d0/11090693500945.d0
      al3N(3,1) = -69155902270455.d0/13319165977508.d0
      al3N(4,1) = 96459639513718.d0/18044420531169.d0
      al3N(5,1) = 680802826223.d0/7936019317924.d0
      al3N(1,2) = -539655284018.d0/6584770896185.d0
      al3N(2,2) = -2297895298714.d0/4249475560743.d0
      al3N(3,2) = -32266676584765.d0/7628952395202.d0
      al3N(4,2) = 47374671169057.d0/8066552623938.d0
      al3N(5,2) = -2052252006800.d0/18582270838531.d0
      al3N(1,3) = 1310161097815.d0/9249264594191.d0
      al3N(2,3) = -16699454119767.d0/22485807752608.d0
      al3N(3,3) = 6500397574583.d0/8403480871568.d0
      al3N(4,3) = 2066487066473.d0/8660385879462.d0
      al3N(5,3) = -1777835010841.d0/17490685973896.d0
      al3N(1,4) = 854410142049.d0/6953488057609.d0
      al3N(2,4) = -2177256661361.d0/6719159179740.d0
      al3N(3,4) = 20357185218447.d0/10240999111918.d0
      al3N(4,4) = -10792423292400.d0/6040590453011.d0
      al3N(1,5) = 417048578066.d0/14947759438935.d0
      al3N(2,5) = 992838867799.d0/14347519519776.d0
      al3N(3,5) = 10653420533222.d0/8563113307551.d0
      al3N(4,5) = -10413366967068.d0/7764183158033.d0
      al3N(1,6) = -117424226259.d0/8877583509224.d0
      al3N(2,6) = 291294039697.d0/4493624958567.d0
      al3N(3,6) = 1950223037569.d0/5991952198682.d0
      al3N(4,6) = -3180846191759.d0/8435679732038.d0
    
    
      al3D(1,1) = 1.d0/1.d0
      al3D(2,1) = 1.d0/1.d0
      al3D(3,1) = 1.d0/1.d0
      al3D(4,1) = 1.d0/1.d0
      al3D(5,1) = 1.d0/1.d0
      al3D(1,2) = 10618085825296.d0/11664012047445.d0
      al3D(2,2) = 10618085825296.d0/11664012047445.d0
      al3D(3,2) = 10618085825296.d0/11664012047445.d0
      al3D(4,2) = 10618085825296.d0/11664012047445.d0
      al3D(5,2) = 10618085825296.d0/11664012047445.d0
      al3D(1,3) = 1222400913289.d0/3949736033075.d0
      al3D(2,3) = 1222400913289.d0/3949736033075.d0
      al3D(3,3) = 1222400913289.d0/3949736033075.d0
      al3D(4,3) = 1222400913289.d0/3949736033075.d0
      al3D(5,3) = 1222400913289.d0/3949736033075.d0
    
    
c   This secondary predictor uses information from the new and old step 
    
    
c   U^{(n+1,4)} = al4_{1}*U^{(n  ,4)} + al4_{2}*U^{(n  ,5)} +
c               + al4_{3}*U^{(n  ,6)} + al4_{4}*U^{(n+1,2)} +
c               + al4_{5}*U^{(n+1,3)}
    
    
c   al4_{i} = \sum_{j=1}^{2*(order)} al4N(i,j)*r^{j} / 
c             \sum_{j=1}^{2*(order)} al4D(i,j)*r^{j} 
    
    
      al4N(1,1) = -1.d0/3120029362199942439971219.d0
      al4N(2,1) = 1.d0/35499000743252678428116.d0
      al4N(3,1) = -61715063944826.d0/7309896013771.d0
      al4N(4,1) = 16937554135729.d0/20578921801152.d0
      al4N(5,1) = 96676337676023.d0/11215847118122.d0
      al4N(1,2) = 50705787112441.d0/14834166081926.d0
      al4N(2,2) = 132087214006161.d0/22678343936063.d0
      al4N(3,2) = -13190803699858.d0/12487528228807.d0
      al4N(4,2) = -23104187274158.d0/7025858887755.d0
      al4N(5,2) = -84858100477173.d0/12806500890004.d0
      al4N(1,3) = -1244726731870.d0/3156207389741.d0
      al4N(2,3) = 383165073545217.d0/11054392434415.d0
      al4N(3,3) = -163507159380887.d0/4984133401290.d0
      al4N(4,3) = -14049590642181.d0/6832537788766.d0
      al4N(5,3) = -8720598832027.d0/9637018600591.d0
      al4N(1,4) = 22797906292193.d0/12665490745416.d0
      al4N(2,4) = 295896270888199.d0/13171456119702.d0
      al4N(3,4) = -142700963964466.d0/6057714649505.d0
      al4N(4,4) = -2394809075965.d0/4953508121693.d0
      al4N(5,4) = -8480134383967.d0/15007887531890.d0
      al4N(1,5) = -1135182128673.d0/3118219248715.d0
      al4N(2,5) = -7898115304907.d0/12732377941718.d0
      al4N(3,5) = 3782054933427.d0/8957305078621.d0
      al4N(4,5) = 2312443347114.d0/8108605632145.d0
      al4N(5,5) = 2734077654837.d0/6237283172479.d0
      al4N(1,6) = 247719513746.d0/9099616704895.d0
      al4N(2,6) = -11338435381549.d0/4738855400152.d0
      al4N(3,6) = 5976315544387.d0/2526524297276.d0
      al4N(1,7) = -861550715461.d0/7595475541482.d0
      al4N(2,7) = -9057135729282.d0/8517138741959.d0
      al4N(3,7) = 18511147585927.d0/15729660898238.d0
    
    
      al4D(1,1) = 1.d0/1.d0
      al4D(2,1) = 1.d0/1.d0
      al4D(3,1) = 1.d0/1.d0
      al4D(4,1) = 1.d0/1.d0
      al4D(5,1) = 1.d0/1.d0
      al4D(1,2) = -15825321138017.d0/9156105750148.d0
      al4D(2,2) = -15825321138017.d0/9156105750148.d0
      al4D(3,2) = -15825321138017.d0/9156105750148.d0
      al4D(4,2) = -15825321138017.d0/9156105750148.d0
      al4D(5,2) = -15825321138017.d0/9156105750148.d0
      al4D(1,3) = -6163363000173.d0/4110832957732.d0
      al4D(2,3) = -6163363000173.d0/4110832957732.d0
      al4D(3,3) = -6163363000173.d0/4110832957732.d0
      al4D(4,3) = -6163363000173.d0/4110832957732.d0
      al4D(5,3) = -6163363000173.d0/4110832957732.d0
      al4D(1,4) = -3604210761703.d0/10587060597062.d0
      al4D(2,4) = -3604210761703.d0/10587060597062.d0
      al4D(3,4) = -3604210761703.d0/10587060597062.d0
      al4D(4,4) = -3604210761703.d0/10587060597062.d0
      al4D(5,4) = -3604210761703.d0/10587060597062.d0
      al4D(1,5) = 2971894851571.d0/18413920481567.d0
      al4D(2,5) = 2971894851571.d0/18413920481567.d0
      al4D(3,5) = 2971894851571.d0/18413920481567.d0
      al4D(4,5) = 2971894851571.d0/18413920481567.d0
      al4D(5,5) = 2971894851571.d0/18413920481567.d0
    
c   Dense output 
    
c   bd_{i} = \sum_{j=1}^{dense_order}  bdense_{ij}*x^{j} 
    
c   U[t^{n} + x*(Delta t)] = U^{(n)}  
c           + (Delta t)\sum_{i=1}^{s} bd_{i}*F^{(i)} 
    
    
      bD(1,1) = 14034291058384.d0/14643615466781.d0
      bD(1,2) = -25215736494927.d0/6674049630709.d0
      bD(1,3) = 41723455951887.d0/9033333068873.d0
      bD(1,4) = -6105373324517.d0/3364520636901.d0
      bD(2,1) = 14034291058384.d0/14643615466781.d0
      bD(2,2) = -25215736494927.d0/6674049630709.d0
      bD(2,3) = 41723455951887.d0/9033333068873.d0
      bD(2,4) = -6105373324517.d0/3364520636901.d0
      bD(3,1) = -203307851416.d0/14003679638037.d0
      bD(3,2) = 26104318035138.d0/6682312545313.d0
      bD(3,3) = -65468352021292.d0/10527533372157.d0
      bD(3,4) = 39018787203283.d0/14374365124454.d0
      bD(4,1) = -18079151498202.d0/13763822068505.d0
      bD(4,2) = 5952586272650.d0/975172309933.d0
      bD(4,3) = -44901778601388.d0/7172115566879.d0
      bD(4,4) = 9433248009269.d0/4784161733256.d0
      bD(5,1) = -14643885851809.d0/8693311047599.d0
      bD(5,2) = 124452067251178.d0/10912632581433.d0
      bD(5,3) = -119066243494329.d0/6546300172030.d0
      bD(5,4) = 91428749833760.d0/10936233599711.d0
      bD(6,1) = 17791496353119.d0/8489256993736.d0
      bD(6,2) = -461495289675690.d0/33300117379213.d0
      bD(6,3) = 225258514035799.d0/10511344076533.d0
      bD(6,4) = -57076810983509.d0/6060951850295.d0
    
      elseif(icase.eq.4) then     !    The scheme you're using is: ESDIRK43I6L2SA
    
      ns = 6

c     The scheme coefficients, A, b, bh, and c  are given by
    
      ai(2,1) = 1.d0/4.d0
      ai(2,2) = 1.d0/4.d0
      ai(3,1) = -1356991263433.d0/26208533697614.d0
      ai(3,2) = -1356991263433.d0/26208533697614.d0
      ai(3,3) = 1.d0/4.d0
      ai(4,1) = -1778551891173.d0/14697912885533.d0
      ai(4,2) = -1778551891173.d0/14697912885533.d0
      ai(4,3) = 7325038566068.d0/12797657924939.d0
      ai(4,4) = 1.d0/4.d0
      ai(5,1) = -24076725932807.d0/39344244018142.d0
      ai(5,2) = -24076725932807.d0/39344244018142.d0
      ai(5,3) = 9344023789330.d0/6876721947151.d0
      ai(5,4) = 11302510524611.d0/18374767399840.d0
      ai(5,5) = 1.d0/4.d0
      ai(6,1) = 657241292721.d0/9909463049845.d0
      ai(6,2) = 657241292721.d0/9909463049845.d0
      ai(6,3) = 1290772910128.d0/5804808736437.d0
      ai(6,4) = 1103522341516.d0/2197678446715.d0
      ai(6,5) = -3.d0/28.d0
      ai(6,6) = 1.d0/4.d0
    
      be(1) = 657241292721.d0/9909463049845.d0
      be(2) = 657241292721.d0/9909463049845.d0
      be(3) = 1290772910128.d0/5804808736437.d0
      be(4) = 1103522341516.d0/2197678446715.d0
      be(5) = -3.d0/28.d0
      be(6) = 1.d0/4.d0
    
      beh(1) = -71925161075.d0/3900939759889.d0
      beh(2) = -71925161075.d0/3900939759889.d0
      beh(3) = 2973346383745.d0/8160025745289.d0
      beh(4) = 3972464885073.d0/7694851252693.d0
      beh(5) = -263368882881.d0/4213126269514.d0
      beh(6) = 3295468053953.d0/15064441987965.d0
    
      bi(1) = 657241292721.d0/9909463049845.d0
      bi(2) = 657241292721.d0/9909463049845.d0
      bi(3) = 1290772910128.d0/5804808736437.d0
      bi(4) = 1103522341516.d0/2197678446715.d0
      bi(5) = -3.d0/28.d0
      bi(6) = 1.d0/4.d0
    
      bih(1) = -71925161075.d0/3900939759889.d0
      bih(2) = -71925161075.d0/3900939759889.d0
      bih(3) = 2973346383745.d0/8160025745289.d0
      bih(4) = 3972464885073.d0/7694851252693.d0
      bih(5) = -263368882881.d0/4213126269514.d0
      bih(6) = 3295468053953.d0/15064441987965.d0
    
c   The primary predictor is now given in two parts; eta and B
    
c   U^{(n+1,i+1)} = \sum_{j=0}^{order}eta_{ij}*r^{j} +  
c                 + \sum_{j=1}^{S-1}\sum_{k=0}^{order}B_{ijk}*r^{k}U^{(n,j+1)} 
    
      svpB(2,6,0) = 1.d0/1.d0
      svpB(3,6,0) = 1.d0/1.d0
      svpB(4,6,0) = 1.d0/1.d0
      svpB(5,6,0) = 1.d0/1.d0
      svpB(6,6,0) = 1.d0/1.d0
      svpB(2,2,1) = 219890023397208985780636766819278.d0/9912605.d0
      svpB(3,2,1) = -122478577601616544434759523644931.d0/26659260.d0
      svpB(4,2,1) = -323433977694424122081221882782924.d0/30122919.d0
      svpB(5,2,1) = -870382228015834449944639591389489.d0/16029356.d0
      svpB(6,2,1) = 255029020026235878970619170641954.d0/43334837.d0
      svpB(2,3,1) = -544185593079871243873493165724806.d0/8261347.d0
      svpB(3,3,1) = 252122423694580958939217258752360.d0/18480809.d0
      svpB(4,3,1) = 909498419497151183375729569865741.d0/28525608.d0
      svpB(5,3,1) = 1112618024927594259053543963498824.d0/6900391.d0
      svpB(6,3,1) = -670190800740870145503910065975448.d0/38350195.d0
      svpB(2,4,1) = 356456114363299350436998535324526.d0/32293867.d0
      svpB(3,4,1) = -107088274061949420821730633814624.d0/46844821.d0
      svpB(4,4,1) = -164665302807579363753716299318970.d0/30820893.d0
      svpB(5,4,1) = -451518457243897268908311545094581.d0/16711417.d0
      svpB(6,4,1) = 132868675319841158650429895784829.d0/45373458.d0
      svpB(2,5,1) = 1287143915553376120291624521065999.d0/122963826.d0
      svpB(3,5,1) = -131435048885155070400711837355553.d0/60627147.d0
      svpB(4,5,1) = -43527590362798483207783456397883.d0/8591006.d0
      svpB(5,5,1) = -576855430836293438668723396859205.d0/22513397.d0
      svpB(6,5,1) = 357796071269309488926250255729507.d0/128840119.d0
      svpB(2,6,1) = -254616070195646261371899055161041.d0/13899463.d0
      svpB(3,6,1) = 350230362807455947735778396088811.d0/92314875.d0
      svpB(4,6,1) = 460438784004149953454409813853129.d0/51929393.d0
      svpB(5,6,1) = 1263957820620523363039101065264769.d0/28188280.d0
      svpB(6,6,1) = -91625255425056985911462181595546.d0/18853531.d0
      svpB(2,2,2) = -198703119498202627540550961047219.d0/31026436.d0
      svpB(3,2,2) = -52924651923786249454041673161825.d0/96331114.d0
      svpB(4,2,2) = -183511855437401026168160018400017.d0/21268539.d0
      svpB(5,2,2) = -776123939217115581532047024906478.d0/30296907.d0
      svpB(6,2,2) = -275459110270418244544663345850953.d0/10752869.d0
      svpB(2,3,2) = 442496873931199756916973296225187.d0/23267995.d0
      svpB(3,3,2) = 131972068369203277573730302589581.d0/80893224.d0
      svpB(4,3,2) = 321868393315018617372569728801229.d0/12562423.d0
      svpB(5,3,2) = 1092882771920344565804882819847709.d0/14366876.d0
      svpB(6,3,2) = 412724309363222487666095882773290.d0/5425613.d0
      svpB(2,4,2) = -296451080665373802212130035815923.d0/93027629.d0
      svpB(3,4,2) = -38416437879510591679265054148019.d0/140526223.d0
      svpB(4,4,2) = -137224503821124270697855526826405.d0/31962218.d0
      svpB(5,4,2) = -217837856962050795076866870102946.d0/17089615.d0
      svpB(6,4,2) = -11747690506515699254076057924938.d0/921619.d0
      svpB(2,5,2) = -559933556124982191117712382147704.d0/185281341.d0
      svpB(3,5,2) = -121338781693785272813237998619799.d0/468032689.d0
      svpB(4,5,2) = -156068395984807348352463679330073.d0/38331559.d0
      svpB(5,5,2) = -696184286544764264359112488890641.d0/57591636.d0
      svpB(6,5,2) = -242527843398508888698703476265394.d0/20063043.d0
      svpB(2,6,2) = 474191161856004125853937489053166.d0/89662449.d0
      svpB(3,6,2) = 112855210041109676765175872992409.d0/248748300.d0
      svpB(4,6,2) = 252265533274139026293295204505187.d0/35404738.d0
      svpB(5,6,2) = 128600261421553955433992987535697.d0/6079096.d0
      svpB(6,6,2) = 220125921626456148040541768298010.d0/10405629.d0
      svpB(2,2,3) = -313650984533464836360476091789580.d0/28532131.d0
      svpB(3,2,3) = 155422408274756960716446161609737.d0/116537999.d0
      svpB(4,2,3) = -992573386149307768151982544819175.d0/85209761.d0
      svpB(5,2,3) = -460441342145171378971976329116813.d0/7853498.d0
      svpB(6,2,3) = -965020045028592507764648559847493.d0/16459823.d0
      svpB(2,3,3) = 248325928545448228682438052488109.d0/7607311.d0
      svpB(3,3,3) = -149090707700153698549296201060589.d0/37646625.d0
      svpB(4,3,3) = 450629043692614305410325806139053.d0/13027691.d0
      svpB(5,3,3) = 2789104564285040253384638337756983.d0/16020467.d0
      svpB(6,3,3) = 2661963752955899720424876940435169.d0/15290177.d0
      svpB(2,4,3) = -157455526827188271921677125760479.d0/28785725.d0
      svpB(3,4,3) = 52136787583464949649823026112411.d0/78565142.d0
      svpB(4,4,3) = -303619254426773027439319405973726.d0/52382683.d0
      svpB(5,4,3) = -3225436978128662968981581825923666.d0/110562817.d0
      svpB(6,4,3) = -495034170011141641034038690311343.d0/16968979.d0
      svpB(2,5,3) = -836340915949157312024117744760259.d0/161227412.d0
      svpB(3,5,3) = 82939056564345049610858798032738.d0/131789557.d0
      svpB(4,5,3) = -300825232498134711450568520713933.d0/54727923.d0
      svpB(5,5,3) = -417636124225492116288742863070117.d0/15095756.d0
      svpB(6,5,3) = -580594337248681343018475417414088.d0/20985997.d0
      svpB(2,6,3) = 492561712601651342589553803682967.d0/54259794.d0
      svpB(3,6,3) = -27062312329249112578080330352067.d0/24572467.d0
      svpB(4,6,3) = 174326492801211141341056827963733.d0/18122581.d0
      svpB(5,6,3) = 286924126962717615048200687904721.d0/5926330.d0
      svpB(6,6,3) = 1188577481745555256797929282885131.d0/24549704.d0
      svpB(2,2,4) = -349358828109776005227648258892020.d0/63371453.d0
      svpB(3,2,4) = 69525710988907584417465643596601.d0/43120751.d0
      svpB(4,2,4) = -125828500524312824955964316572381.d0/45289844.d0
      svpB(5,2,4) = -378717519916246633736000895039860.d0/11472271.d0
      svpB(6,2,4) = -206820767349578528742276622490729.d0/6265102.d0
      svpB(2,3,4) = 213203694830563685380374718407537.d0/13023815.d0
      svpB(3,3,4) = -187893348068300143487123081063737.d0/39244051.d0
      svpB(4,3,4) = 69394299743624823173033375758837.d0/8411378.d0
      svpB(5,3,4) = 514253245449645083771450933115185.d0/5246052.d0
      svpB(6,3,4) = 2190688582413307421885781145402282.d0/22347873.d0
      svpB(2,4,4) = -69628694874163149555848202855252.d0/25382939.d0
      svpB(3,4,4) = 165095628005822425326061986651324.d0/205782367.d0
      svpB(4,4,4) = -222439312659942559458619630287439.d0/160903339.d0
      svpB(5,4,4) = -469441079787050489388403540036316.d0/28578985.d0
      svpB(6,4,4) = -5724895745308782726503069386371106.d0/348524483.d0
      svpB(2,5,4) = -122328480717950167496025601013164.d0/47023777.d0
      svpB(3,5,4) = 21189751600757765537913298413575.d0/27850611.d0
      svpB(4,5,4) = -101646488730627532268938033975634.d0/77532211.d0
      svpB(5,5,4) = -595928220326562851771870398135260.d0/38255681.d0
      svpB(6,5,4) = -254196831714793321676053207565289.d0/16318195.d0
      svpB(2,6,4) = 186674929859490704523474154745793.d0/41005100.d0
      svpB(3,6,4) = -80717871041767679958736397110633.d0/60623431.d0
      svpB(4,6,4) = 98263807293398089889122359181234.d0/42829729.d0
      svpB(5,6,4) = 312188369926696791141431411910165.d0/11451982.d0
      svpB(6,6,4) = 648585170435647738080354656718649.d0/23792000.d0
    
c   U^{(n+1,i+1)} = \sum_{j=0}^{order}eta_{ij}*r^{j} +  
c                 + \sum_{j=1}^{S-1}\sum_{k=0}^{order}B_{ijk}*r^{k}U^{(n,j+1)} 
    
      svpB(2,1,1) = 422745031125402318328198278919132.d0/10437819.d0
      svpB(3,1,1) = -376723018524318252955946898826837.d0/44911658.d0
      svpB(4,1,1) = -486581251188215455505512476362524.d0/24820799.d0
      svpB(5,1,1) = -1026975577001199116417514137993217.d0/10358936.d0
      svpB(6,1,1) = 465209695446449077604652998879822.d0/43295753.d0
      svpB(2,1,2) = -311870576472470823194496561915413.d0/26671693.d0
      svpB(3,1,2) = -86789196956017070034461460007243.d0/86521336.d0
      svpB(4,1,2) = -292753744229461686687460730416054.d0/18583391.d0
      svpB(5,1,2) = -1988643965481905762329478852758123.d0/42518039.d0
      svpB(6,1,2) = -675852271921322622458764498797685.d0/14450004.d0
      svpB(2,1,3) = -1047851454510196851967869948989606.d0/52207903.d0
      svpB(3,1,3) = 79439518040024020446321683871439.d0/32624172.d0
      svpB(4,1,3) = -258053628050931174835416973311343.d0/12133488.d0
      svpB(5,1,3) = -2418086410811277148670982492013041.d0/22589656.d0
      svpB(6,1,3) = -7740474878463104185760691786782173.d0/72311173.d0
      svpB(2,1,4) = -184946735025212855766530192201418.d0/18374587.d0
      svpB(3,1,4) = 166140109103758108262055638581542.d0/56437055.d0
      svpB(4,1,4) = -139690351960606895964677865908967.d0/27538304.d0
      svpB(5,1,4) = -317667991019401582935160785235891.d0/5270556.d0
      svpB(6,1,4) = -857223134013407616832609712266217.d0/14222530.d0
    
c  The secondary predictor is used well into the next step
c  It uses only information from the new step 
    
c   Z_{0}^{(n+1,i)} = \sum_{j=2}^{i-1} \alpha_{ij} Z_{0}^{(n+1,j)} 
    
      alpha(3,2) = 1513744654945.d0/5168247530883.d0
      alpha(4,2) = 69831889866857.d0/49020727896598.d0
      alpha(4,3) = -10020181211492.d0/11124427962077.d0
      alpha(5,2) = 32541972357828.d0/11516941232453.d0
      alpha(5,3) = -72846998129358.d0/8682134072279.d0
      alpha(5,4) = 13173520290768.d0/9369724731565.d0
      alpha(6,2) = 4849389959219.d0/4004581835792.d0
      alpha(6,3) = -52662912326633.d0/14645239750818.d0
      alpha(6,4) = 5012788920243.d0/8319192819736.d0
      alpha(6,5) = 320790218736368104353345.d0/561382882788644182618354.d0
    
c   This secondary predictor uses information from the new and old step 
    
    
c   U^{(n+1,3)} = al3_{1}*U^{(n  ,3)} + al3_{2}*U^{(n  ,4)} +
c               + al3_{3}*U^{(n  ,5)} + al3_{4}*U^{(n  ,6)} +
c               + al3_{5}*U^{(n+1,2)}
    
    
c   al3_{i} = \sum_{j=1}^{2*(order)} al3N(i,j)*r^{j} / 
c             \sum_{j=1}^{2*(order)} al3D(i,j)*r^{j} 
    
    
      al3N(1,1) = -1.d0/10890451872054993415699875.d0
      al3N(2,1) = 1.d0/2632328084699137051343892.d0
      al3N(3,1) = -9581044292864.d0/3036130422891.d0
      al3N(4,1) = 37371565715578.d0/9182451754533.d0
      al3N(5,1) = 156753391512.d0/1827251437969.d0
      al3N(1,2) = 517821742308.d0/14960342112319.d0
      al3N(2,2) = -3304844250148.d0/10419794940325.d0
      al3N(3,2) = -27907988989060.d0/10040563442779.d0
      al3N(4,2) = 30859595811274.d0/8385462076605.d0
      al3N(5,2) = -920787896312.d0/13769949346089.d0
      al3N(1,3) = 1451121336842.d0/11219433578639.d0
      al3N(2,3) = -9431042007137.d0/17624670256220.d0
      al3N(3,3) = 2535600533992.d0/4969611560057.d0
      al3N(4,3) = -9225562886641.d0/88319696830030.d0
      al3N(5,3) = -1.d0/24994550938696907515378839.d0
      al3N(1,4) = 490942363772.d0/6371408565119.d0
      al3N(2,4) = -14884122899666.d0/94967553852255.d0
      al3N(3,4) = 10056373166971.d0/7908786826820.d0
      al3N(4,4) = -11344222220854.d0/9518004467609.d0
      al3N(1,5) = 1.d0/29359789295255660894023437.d0
      al3N(2,5) = -1.d0/14434441778013562049492968.d0
      al3N(3,5) = 6529853183329.d0/7882230057749.d0
      al3N(4,5) = -11545961711713.d0/13937208677529.d0
      al3N(3,6) = 765082525703.d0/3351122783842.d0
      al3N(4,6) = -765082525703.d0/3351122783842.d0
    
    
      al3D(1,1) = 1.d0/1.d0
      al3D(2,1) = 1.d0/1.d0
      al3D(3,1) = 1.d0/1.d0
      al3D(4,1) = 1.d0/1.d0
      al3D(5,1) = 1.d0/1.d0
      al3D(1,2) = 9311950492833.d0/16894565644852.d0
      al3D(2,2) = 9311950492833.d0/16894565644852.d0
      al3D(3,2) = 9311950492833.d0/16894565644852.d0
      al3D(4,2) = 9311950492833.d0/16894565644852.d0
      al3D(5,2) = 9311950492833.d0/16894565644852.d0
      al3D(1,3) = 1.d0/8208888499118412047468694.d0
      al3D(2,3) = 1.d0/8208888499118412047468694.d0
      al3D(3,3) = 1.d0/8208888499118412047468694.d0
      al3D(4,3) = 1.d0/8208888499118412047468694.d0
      al3D(5,3) = 1.d0/8208888499118412047468694.d0
    
    
c   This secondary predictor uses information from the new and old step 
    
    
c   U^{(n+1,4)} = al4_{1}*U^{(n  ,4)} + al4_{2}*U^{(n  ,5)} +
c               + al4_{3}*U^{(n  ,6)} + al4_{4}*U^{(n+1,2)} +
c               + al4_{5}*U^{(n+1,3)}
    
    
c   al4_{i} = \sum_{j=1}^{2*(order)} al4N(i,j)*r^{j} / 
c             \sum_{j=1}^{2*(order)} al4D(i,j)*r^{j} 
    
    
      al4N(2,1) = -1.d0/19056409840493564021771.d0
      al4N(3,1) = 900742089531.d0/1891525494161.d0
      al4N(4,1) = 8824481787402.d0/6194627144611.d0
      al4N(5,1) = -4162625043215.d0/4621355786787.d0
      al4N(1,2) = -1.d0/162303052261891443957616.d0
      al4N(2,2) = -122729415331101.d0/5491766346854.d0
      al4N(3,2) = 403969549264835.d0/19143128798259.d0
      al4N(4,2) = 116769401198809.d0/28688488357637.d0
      al4N(5,2) = 6597687232683.d0/7235859154115.d0
      al4N(1,3) = -1.d0/2585939146765945889065156.d0
      al4N(2,3) = -2873832040757861.d0/47349169357959.d0
      al4N(3,3) = 375545749522427.d0/6577558605251.d0
      al4N(4,3) = 25274804181918.d0/9244575625145.d0
      al4N(5,3) = 8474694016493.d0/2741205273822.d0
      al4N(1,4) = 21460186801460.d0/17062836057241.d0
      al4N(2,4) = -384616553535995.d0/11195183959106.d0
      al4N(3,4) = 273105493418396.d0/8251464302231.d0
      al4N(4,4) = 1.d0/899376204351444295686985.d0
      al4N(5,4) = 1.d0/401116190282873601902646.d0
      al4N(1,5) = 1.d0/1199152589840250996714968.d0
      al4N(2,5) = 35703418944685.d0/12549453889854.d0
      al4N(3,5) = -16591003498784.d0/5831599313135.d0
      al4N(2,6) = 51820114543753.d0/9748759262877.d0
      al4N(3,6) = -51820114543753.d0/9748759262877.d0
      al4N(2,7) = 22176242702234.d0/10504320017751.d0
      al4N(3,7) = -22176242702234.d0/10504320017751.d0
    
    
      al4D(1,1) = 1.d0/1.d0
      al4D(2,1) = 1.d0/1.d0
      al4D(3,1) = 1.d0/1.d0
      al4D(4,1) = 1.d0/1.d0
      al4D(5,1) = 1.d0/1.d0
      al4D(1,2) = 57781619523892.d0/15463068177921.d0
      al4D(2,2) = 57781619523892.d0/15463068177921.d0
      al4D(3,2) = 57781619523892.d0/15463068177921.d0
      al4D(4,2) = 57781619523892.d0/15463068177921.d0
      al4D(5,2) = 57781619523892.d0/15463068177921.d0
      al4D(1,3) = 56356324838461.d0/25315491639988.d0
      al4D(2,3) = 56356324838461.d0/25315491639988.d0
      al4D(3,3) = 56356324838461.d0/25315491639988.d0
      al4D(4,3) = 56356324838461.d0/25315491639988.d0
      al4D(5,3) = 56356324838461.d0/25315491639988.d0
      al4D(1,4) = 1.d0/1016229250141687978508712.d0
      al4D(2,4) = 1.d0/1016229250141687978508712.d0
      al4D(3,4) = 1.d0/1016229250141687978508712.d0
      al4D(4,4) = 1.d0/1016229250141687978508712.d0
      al4D(5,4) = 1.d0/1016229250141687978508712.d0
    
c   Dense output 
    
c   bd_{i} = \sum_{j=1}^{dense_order}  bdense_{ij}*x^{j} 
    
c   U[t^{n} + x*(Delta t)] = U^{(n)}  
c           + (Delta t)\sum_{i=1}^{s} bd_{i}*F^{(i)} 
    
    
      bD(1,1) = 13296219569665.d0/11073869861459.d0
      bD(1,2) = -46101377119439.d0/8039770269915.d0
      bD(1,3) = 55943395092913.d0/6879771013856.d0
      bD(1,4) = -34619626224326.d0/9802340079121.d0
      bD(2,1) = 13296219569665.d0/11073869861459.d0
      bD(2,2) = -46101377119439.d0/8039770269915.d0
      bD(2,3) = 55943395092913.d0/6879771013856.d0
      bD(2,4) = -34619626224326.d0/9802340079121.d0
      bD(3,1) = -2276017419075.d0/5117804105921.d0
      bD(3,2) = 87682662658264.d0/11731837202815.d0
      bD(3,3) = -239878786873189.d0/18852188973841.d0
      bD(3,4) = 73302989622937.d0/12387768562097.d0
      bD(4,1) = -19596115037837.d0/19512953221608.d0
      bD(4,2) = 84802285654569.d0/20667355910242.d0
      bD(4,3) = -10997452681547.d0/3452792107192.d0
      bD(4,4) = 1852267820943.d0/3148600010176.d0
      bD(5,1) = -81442679768748405998409084673220.d0/20250073.d0
      bD(5,2) = 89418457410857731755734170268813.d0/3439292.d0
      bD(5,3) = -1227712638136562287348854628651585.d0/30744591.d0
      bD(5,4) = 2310540504707258658998512974317872.d0/128682235.d0
      bD(6,1) = 455859306115309670199370074739407.d0/113345782.d0
      bD(6,2) = -1031800557941607588560248870026667.d0/39686028.d0
      bD(6,3) = 173869383502071871725018731881568.d0/4354067.d0
      bD(6,4) = -790819010291115811425734981828118.d0/44043529.d0
    
      elseif(icase.eq.5) then     !    The scheme you're using is: ESDIRK536L2SA

      ns = 6
    
c     The scheme coefficients, A, b, bh, and c  are given by
    
      ai(2,1) = 3282482714977.d0/11805205429139.d0
      ai(2,2) = 3282482714977.d0/11805205429139.d0
      ai(3,1) = 606638434273.d0/1934588254988.d0
      ai(3,2) = 2719561380667.d0/6223645057524.d0
      ai(3,3) = 3282482714977.d0/11805205429139.d0
      ai(4,1) = -651839358321.d0/6893317340882.d0
      ai(4,2) = -1510159624805.d0/11312503783159.d0
      ai(4,3) = 235043282255.d0/4700683032009.d0
      ai(4,4) = 3282482714977.d0/11805205429139.d0
      ai(5,1) = -5266892529762.d0/23715740857879.d0
      ai(5,2) = -1007523679375.d0/10375683364751.d0
      ai(5,3) = 521543607658.d0/16698046240053.d0
      ai(5,4) = 514935039541.d0/7366641897523.d0
      ai(5,5) = 3282482714977.d0/11805205429139.d0
      ai(6,1) = -6225479754948.d0/6925873918471.d0
      ai(6,2) = 6894665360202.d0/11185215031699.d0
      ai(6,3) = -2508324082331.d0/20512393166649.d0
      ai(6,4) = -7289596211309.d0/4653106810017.d0
      ai(6,5) = 39811658682819.d0/14781729060964.d0
      ai(6,6) = 3282482714977.d0/11805205429139.d0
    
      be(1) = -6225479754948.d0/6925873918471.d0
      be(2) = 6894665360202.d0/11185215031699.d0
      be(3) = -2508324082331.d0/20512393166649.d0
      be(4) = -7289596211309.d0/4653106810017.d0
      be(5) = 39811658682819.d0/14781729060964.d0
      be(6) = 3282482714977.d0/11805205429139.d0
    
      beh(1) = -2512930284403.d0/5616797563683.d0
      beh(2) = 5849584892053.d0/8244045029872.d0
      beh(3) = -718651703996.d0/6000050726475.d0
      beh(4) = -18982822128277.d0/13735826808854.d0
      beh(5) = 23127941173280.d0/11608435116569.d0
      beh(6) = 2847520232427.d0/11515777524847.d0
    
      bi(1) = -6225479754948.d0/6925873918471.d0
      bi(2) = 6894665360202.d0/11185215031699.d0
      bi(3) = -2508324082331.d0/20512393166649.d0
      bi(4) = -7289596211309.d0/4653106810017.d0
      bi(5) = 39811658682819.d0/14781729060964.d0
      bi(6) = 3282482714977.d0/11805205429139.d0
    
      bih(1) = -2512930284403.d0/5616797563683.d0
      bih(2) = 5849584892053.d0/8244045029872.d0
      bih(3) = -718651703996.d0/6000050726475.d0
      bih(4) = -18982822128277.d0/13735826808854.d0
      bih(5) = 23127941173280.d0/11608435116569.d0
      bih(6) = 2847520232427.d0/11515777524847.d0
    
c   The primary predictor is now given in two parts; eta and B
    
c   U^{(n+1,i+1)} = \sum_{j=0}^{order}eta_{ij}*r^{j} +  
c                 + \sum_{j=1}^{S-1}\sum_{k=0}^{order}B_{ijk}*r^{k}U^{(n,j+1)} 
    
      svpB(2,6,0) = 1.d0/1.d0
      svpB(3,6,0) = 1.d0/1.d0
      svpB(4,6,0) = 1.d0/1.d0
      svpB(5,6,0) = 1.d0/1.d0
      svpB(6,6,0) = 1.d0/1.d0
      svpB(2,2,1) = -82542913125771.d0/11431673317210.d0
      svpB(3,2,1) = -90481596039380.d0/8371413096451.d0
      svpB(4,2,1) = 32656904234237.d0/60932980162682.d0
      svpB(5,2,1) = 11101741581641.d0/4588087965044.d0
      svpB(6,2,1) = 52055568304309.d0/10921319491275.d0
      svpB(2,3,1) = -20158895736830.d0/5667567013989.d0
      svpB(3,3,1) = -24490327514075.d0/6278848881427.d0
      svpB(4,3,1) = 20276826171115.d0/15725178722058.d0
      svpB(5,3,1) = 52715726415351.d0/17689224959852.d0
      svpB(6,3,1) = 253547522099423.d0/20662642022034.d0
      svpB(2,4,1) = 332634389478349.d0/6857464514739.d0
      svpB(3,4,1) = 715681272169049.d0/10785816721495.d0
      svpB(4,4,1) = -102465991997933.d0/12640981948807.d0
      svpB(5,4,1) = -30841448330786.d0/1279106068349.d0
      svpB(6,4,1) = -434094266623475.d0/5740635680339.d0
      svpB(2,5,1) = -1325806042277722.d0/23050638706983.d0
      svpB(3,5,1) = -902187258931513.d0/11441847638077.d0
      svpB(4,5,1) = 90679803178185.d0/9556621483013.d0
      svpB(5,5,1) = 194871931011478.d0/6867421279315.d0
      svpB(6,5,1) = 4362554610984755.d0/49308216454132.d0
      svpB(2,6,1) = 34024598571019.d0/4981302304128.d0
      svpB(3,6,1) = 97841758302109.d0/10696785346317.d0
      svpB(4,6,1) = -6492155192027.d0/5059703252592.d0
      svpB(5,6,1) = -27238487085053.d0/7478271113599.d0
      svpB(6,6,1) = -73874347288060.d0/6146385939041.d0
      svpB(2,2,2) = -43515166245932.d0/10156364104013.d0
      svpB(3,2,2) = -119323244159987.d0/8140411642060.d0
      svpB(4,2,2) = -918395986211.d0/6628959734825.d0
      svpB(5,2,2) = -489301283815.d0/9810457451262.d0
      svpB(6,2,2) = -128171615698780.d0/9251395828939.d0
      svpB(2,3,2) = -668941988719.d0/12127821754347.d0
      svpB(3,3,2) = -1347747827837.d0/7142112527480.d0
      svpB(4,3,2) = -22635797931.d0/12691353648070.d0
      svpB(5,3,2) = -5818326431.d0/9061659103085.d0
      svpB(6,3,2) = -1670050926663.d0/9363578427929.d0
      svpB(2,4,2) = 143694143079734.d0/8652986972701.d0
      svpB(3,4,2) = 540670632502610.d0/9516636106659.d0
      svpB(4,4,2) = 4619126697991.d0/8602103169535.d0
      svpB(5,4,2) = 1445799248285.d0/7479115948712.d0
      svpB(6,4,2) = 376196027145406.d0/7005820037111.d0
      svpB(2,5,2) = -111452950884418.d0/8405881167591.d0
      svpB(3,5,2) = -360841182434856.d0/7954840411631.d0
      svpB(4,5,2) = -5310684525521.d0/12386824011787.d0
      svpB(5,5,2) = -3378858926069.d0/21891575724342.d0
      svpB(6,5,2) = -259091081405354.d0/6043129869547.d0
      svpB(2,6,2) = 36812628397003.d0/23383551672302.d0
      svpB(3,6,2) = 73589781094229.d0/13663276571870.d0
      svpB(4,6,2) = 1078711484166.d0/21190282202519.d0
      svpB(5,6,2) = 307368774457.d0/16772147878360.d0
      svpB(6,6,2) = 47671944612821.d0/9364709417825.d0
      svpB(2,2,3) = -6864976103930.d0/2820239732711.d0
      svpB(3,2,3) = -196685983120715.d0/16183894195902.d0
      svpB(4,2,3) = -3523541908725.d0/8644154229086.d0
      svpB(5,2,3) = -1683853253387.d0/12613312762304.d0
      svpB(6,2,3) = -82614555427873.d0/8755309641733.d0
      svpB(2,3,3) = 735341690159.d0/9088178865562.d0
      svpB(3,3,3) = 4628134577768.d0/11456613935125.d0
      svpB(4,3,3) = 183087247429.d0/13512690854148.d0
      svpB(5,3,3) = 41982510619.d0/9460933955384.d0
      svpB(6,3,3) = 3742961739566.d0/11933587109139.d0
      svpB(2,4,3) = 1026805804406.d0/9875325053413.d0
      svpB(3,4,3) = 166645147640.d0/321009865539.d0
      svpB(4,4,3) = 239463621650.d0/13753039799837.d0
      svpB(5,4,3) = 38518200477.d0/6754714558021.d0
      svpB(6,4,3) = 3147188510937.d0/7808251345310.d0
      svpB(2,5,3) = 43864831954784.d0/4165610514137.d0
      svpB(3,5,3) = 358424395702205.d0/6817456495839.d0
      svpB(4,5,3) = 31938465645796.d0/18112245993865.d0
      svpB(5,5,3) = 5723374838799.d0/9910425661999.d0
      svpB(6,5,3) = 731674362045510.d0/17924538030487.d0
      svpB(2,6,3) = 1534971613487.d0/2443328683625.d0
      svpB(3,6,3) = 23937861909488.d0/7631844711943.d0
      svpB(4,6,3) = 2298433266851.d0/21847890606660.d0
      svpB(5,6,3) = 234264708324.d0/6799333448275.d0
      svpB(6,6,3) = 15227111379205.d0/6252691654241.d0
      svpB(2,2,4) = -1516947966263.d0/12489081406224.d0
      svpB(3,2,4) = -32517594340658.d0/14951915149105.d0
      svpB(4,2,4) = -461679922521.d0/3036688493927.d0
      svpB(5,2,4) = -1087752780311.d0/11449212403488.d0
      svpB(6,2,4) = -6179379518907.d0/2782070688445.d0
      svpB(2,3,4) = 5110600098529.d0/10967594368181.d0
      svpB(3,3,4) = 25911835297317.d0/10783594984817.d0
      svpB(4,3,4) = 2861494983756.d0/12221565110893.d0
      svpB(5,3,4) = 493780100645.d0/4242766019489.d0
      svpB(6,3,4) = 5490496843504.d0/7148411431219.d0
      svpB(2,4,4) = -133145523865610.d0/17928530836257.d0
      svpB(3,4,4) = -532278564674763.d0/11202091419503.d0
      svpB(4,4,4) = -4890448291220.d0/1144397320311.d0
      svpB(5,4,4) = -33551204507569.d0/14978927362272.d0
      svpB(6,4,4) = -179095810998335.d0/7378057762217.d0
      svpB(2,5,4) = 195846562818182.d0/15512095225703.d0
      svpB(3,5,4) = 1453082124783817.d0/16780106196792.d0
      svpB(4,5,4) = 103905394096777.d0/13659574574970.d0
      svpB(5,5,4) = 55431625820011.d0/13684011553818.d0
      svpB(6,5,4) = 99695054366603.d0/2040529347643.d0
      svpB(2,6,4) = -3519775098398.d0/8250127743109.d0
      svpB(3,6,4) = -9254248789401.d0/5423506364530.d0
      svpB(4,6,4) = -4348287350321.d0/23459527968234.d0
      svpB(5,6,4) = -703662240977.d0/8188133033045.d0
      svpB(6,6,4) = -782851576921.d0/13295627114126.d0
    
c   U^{(n+1,i+1)} = \sum_{j=0}^{order}eta_{ij}*r^{j} +  
c                 + \sum_{j=1}^{S-1}\sum_{k=0}^{order}B_{ijk}*r^{k}U^{(n,j+1)} 
    
      svpB(2,1,1) = 149962392313113.d0/11573699697721.d0
      svpB(3,1,1) = 249934669427311.d0/13840756915724.d0
      svpB(4,1,1) = -34932867438334.d0/18145821647913.d0
      svpB(5,1,1) = -93398288179681.d0/15509489304649.d0
      svpB(6,1,1) = -207246017712989.d0/11593868174950.d0
      svpB(2,1,2) = -4733166553849.d0/8132643804355.d0
      svpB(3,1,2) = -6297482223884.d0/3162791872587.d0
      svpB(4,1,2) = -251063352432.d0/13340775542591.d0
      svpB(5,1,2) = -56324151705.d0/8313615072368.d0
      svpB(6,1,2) = -16826173193803.d0/8940938517902.d0
      svpB(2,1,3) = -49091255908977.d0/5510196836333.d0
      svpB(3,1,3) = -527609952128295.d0/11861471556641.d0
      svpB(4,1,3) = -15707952923892.d0/10528794031589.d0
      svpB(5,1,3) = -5043846409658.d0/10322924500043.d0
      svpB(6,1,3) = -683218558963232.d0/19782927898711.d0
      svpB(2,1,4) = -57133750170775.d0/11165853871844.d0
      svpB(3,1,4) = -501318107023375.d0/13332471783952.d0
      svpB(4,1,4) = -25752460527805.d0/7972537389761.d0
      svpB(5,1,4) = -6232842630574.d0/3569015432749.d0
      svpB(6,1,4) = -256138048589312.d0/11101981418429.d0
    
c  The secondary predictor is used well into the next step
c  It uses only information from the new step 
    
c   Z_{0}^{(n+1,i)} = \sum_{j=2}^{i-1} \alpha_{ij} Z_{0}^{(n+1,j)} 
    
      alpha(3,2) = 18421317152934.d0/9959385188353.d0
      alpha(4,2) = 4964164777445.d0/14046616225789.d0
      alpha(4,3) = -675475855640.d0/7197561463617.d0
      alpha(5,2) = 357000086811.d0/7630286006497.d0
      alpha(5,3) = -185173654406.d0/12396873177479.d0
      alpha(5,4) = 9337281144505.d0/18922217055058.d0
      alpha(6,2) = 14240387063408.d0/11311716834071.d0
      alpha(6,3) = 4579748402515.d0/7640495043718.d0
      alpha(6,4) = -11027417373738.d0/9784332577681.d0
      alpha(6,5) = -44912392519429.d0/13214010916750.d0
    
c   This secondary predictor uses information from the new and old step 
    
    
c   U^{(n+1,3)} = al3_{1}*U^{(n  ,3)} + al3_{2}*U^{(n  ,4)} +
c               + al3_{3}*U^{(n  ,5)} + al3_{4}*U^{(n  ,6)} +
c               + al3_{5}*U^{(n+1,2)}
    
    
c   al3_{i} = \sum_{j=1}^{2*(order)} al3N(i,j)*r^{j} / 
c             \sum_{j=1}^{2*(order)} al3D(i,j)*r^{j} 
    
    
      al3N(1,1) = 89245082273845.d0/22014966114964.d0
      al3N(2,1) = -96324228826719.d0/7328376125912.d0
      al3N(3,1) = 135366616472117.d0/11238019090519.d0
      al3N(4,1) = -55169978766341.d0/10261454413089.d0
      al3N(5,1) = 13935465320841.d0/4073288530492.d0
      al3N(1,2) = 88615630371983.d0/6930676755699.d0
      al3N(2,2) = -581823945545452.d0/37473907056749.d0
      al3N(3,2) = 180864699013328.d0/11175394400197.d0
      al3N(4,2) = -159255696351983.d0/12408668498444.d0
      al3N(5,2) = -1452573603951.d0/1904784819457.d0
      al3N(1,3) = 33489083448362.d0/2736067611895.d0
      al3N(2,3) = 20080266694469.d0/8202745577789.d0
      al3N(3,3) = -25917409584169.d0/11491178746622.d0
      al3N(4,3) = -78815242482718.d0/7722919213869.d0
      al3N(5,3) = -32494235400066.d0/12367870112575.d0
      al3N(1,4) = 77730424470311.d0/8398958099016.d0
      al3N(2,4) = -87257518317950.d0/3905793589657.d0
      al3N(3,4) = 1140963759355157.d0/55507546436016.d0
      al3N(4,4) = -56568761553235.d0/7573456213892.d0
      al3N(1,5) = 49340257307215.d0/11908486562732.d0
      al3N(2,5) = -53366826640609.d0/7005013133591.d0
      al3N(3,5) = 101611529327682.d0/13693801343939.d0
      al3N(4,5) = -61677220046698.d0/15633617644133.d0
      al3N(1,6) = 15712997340905.d0/24019903971643.d0
      al3N(2,6) = 6720016475445.d0/13349196787414.d0
      al3N(3,6) = -29046261815468.d0/62860256764563.d0
      al3N(4,6) = -3595517636558.d0/5169752278891.d0
    
    
      al3D(1,1) = 1.d0/1.d0
      al3D(2,1) = 1.d0/1.d0
      al3D(3,1) = 1.d0/1.d0
      al3D(4,1) = 1.d0/1.d0
      al3D(5,1) = 1.d0/1.d0
      al3D(1,2) = -1211456624021.d0/7931457830749.d0
      al3D(2,2) = -1211456624021.d0/7931457830749.d0
      al3D(3,2) = -1211456624021.d0/7931457830749.d0
      al3D(4,2) = -1211456624021.d0/7931457830749.d0
      al3D(5,2) = -1211456624021.d0/7931457830749.d0
      al3D(1,3) = -3864710084597.d0/9655868565109.d0
      al3D(2,3) = -3864710084597.d0/9655868565109.d0
      al3D(3,3) = -3864710084597.d0/9655868565109.d0
      al3D(4,3) = -3864710084597.d0/9655868565109.d0
      al3D(5,3) = -3864710084597.d0/9655868565109.d0
    
    
c   This secondary predictor uses information from the new and old step 
    
    
c   U^{(n+1,4)} = al4_{1}*U^{(n  ,4)} + al4_{2}*U^{(n  ,5)} +
c               + al4_{3}*U^{(n  ,6)} + al4_{4}*U^{(n+1,2)} +
c               + al4_{5}*U^{(n+1,3)}
    
    
c   al4_{i} = \sum_{j=1}^{2*(order)} al4N(i,j)*r^{j} / 
c             \sum_{j=1}^{2*(order)} al4D(i,j)*r^{j} 
    
    
      al4N(1,1) = -1.d0/7821778919366095621289731.d0
      al4N(2,1) = 1.d0/8532498584138126038236551.d0
      al4N(3,1) = 27320465403698.d0/20870012755029.d0
      al4N(4,1) = -4178400388143.d0/9283508650649.d0
      al4N(5,1) = 1008948635318.d0/7155109516093.d0
      al4N(1,2) = -12527984128842.d0/6681296856071.d0
      al4N(2,2) = 1345103814463.d0/852985133481.d0
      al4N(3,2) = 8613927114691.d0/6411592064617.d0
      al4N(4,2) = -1172761859177.d0/7804601630820.d0
      al4N(5,2) = 1388949887156.d0/22328757657645.d0
      al4N(1,3) = -40395255197581.d0/13416612375313.d0
      al4N(2,3) = 31840062916700.d0/11536049795257.d0
      al4N(3,3) = 2393484461908.d0/2516998101195.d0
      al4N(4,3) = -586567855645.d0/15214843668582.d0
      al4N(5,3) = 519553947843.d0/8348342517797.d0
      al4N(1,4) = -4276117350401.d0/8435924806167.d0
      al4N(2,4) = 2363899077888.d0/5211253904195.d0
      al4N(3,4) = 5036285429629.d0/9975741656383.d0
      al4N(4,4) = -124545027218.d0/674954668705.d0
      al4N(5,4) = 788948539998.d0/13841596270987.d0
      al4N(1,5) = -9516961543167.d0/13637204283686.d0
      al4N(2,5) = 2879423332205.d0/4906125295378.d0
      al4N(3,5) = 3385231986717.d0/16286290056250.d0
      al4N(4,5) = -407733697507.d0/6872203511694.d0
      al4N(5,5) = 166356948889.d0/12232598369860.d0
      al4N(1,6) = -7037490534495.d0/11173162418723.d0
      al4N(2,6) = 6616883883930.d0/11459943858107.d0
      al4N(3,6) = 765190151792.d0/14584970430439.d0
      al4N(1,7) = -709175884961.d0/5709557633834.d0
      al4N(2,7) = 1338946990865.d0/11258932867797.d0
      al4N(3,7) = 99659436110.d0/18855358307951.d0
    
    
      al4D(1,1) = 1.d0/1.d0
      al4D(2,1) = 1.d0/1.d0
      al4D(3,1) = 1.d0/1.d0
      al4D(4,1) = 1.d0/1.d0
      al4D(5,1) = 1.d0/1.d0
      al4D(1,2) = 6262327986509.d0/6541751770438.d0
      al4D(2,2) = 6262327986509.d0/6541751770438.d0
      al4D(3,2) = 6262327986509.d0/6541751770438.d0
      al4D(4,2) = 6262327986509.d0/6541751770438.d0
      al4D(5,2) = 6262327986509.d0/6541751770438.d0
      al4D(1,3) = 13298782186927.d0/18373036405755.d0
      al4D(2,3) = 13298782186927.d0/18373036405755.d0
      al4D(3,3) = 13298782186927.d0/18373036405755.d0
      al4D(4,3) = 13298782186927.d0/18373036405755.d0
      al4D(5,3) = 13298782186927.d0/18373036405755.d0
      al4D(1,4) = 1505067413721.d0/4644573084725.d0
      al4D(2,4) = 1505067413721.d0/4644573084725.d0
      al4D(3,4) = 1505067413721.d0/4644573084725.d0
      al4D(4,4) = 1505067413721.d0/4644573084725.d0
      al4D(5,4) = 1505067413721.d0/4644573084725.d0
      al4D(1,5) = 300505471121.d0/5873531474952.d0
      al4D(2,5) = 300505471121.d0/5873531474952.d0
      al4D(3,5) = 300505471121.d0/5873531474952.d0
      al4D(4,5) = 300505471121.d0/5873531474952.d0
      al4D(5,5) = 300505471121.d0/5873531474952.d0
    
c   Dense output 
    
c   bd_{i} = \sum_{j=1}^{dense_order}  bdense_{ij}*x^{j} 
    
c   U[t^{n} + x*(Delta t)] = U^{(n)}  
c           + (Delta t)\sum_{i=1}^{s} bd_{i}*F^{(i)} 
    
    
      bD(1,1) = 21792830975905.d0/10173202522378.d0
      bD(1,2) = -98958054011471.d0/6239518801079.d0
      bD(1,3) = 268209432355561.d0/12361169151423.d0
      bD(1,4) = -85334073676666.d0/9610875823743.d0
      bD(2,1) = -1139554570237.d0/9034749720566.d0
      bD(2,2) = 4477284553623.d0/16023590656072.d0
      bD(2,3) = 25513075105760.d0/11164534871043.d0
      bD(2,4) = -14314284012103.d0/7856062287211.d0
      bD(3,1) = -2131074325779.d0/9426370871165.d0
      bD(3,2) = 7592368030783.d0/12279056716840.d0
      bD(3,3) = -8205496983035.d0/7833088415518.d0
      bD(3,4) = 5317037053385.d0/9975362047644.d0
      bD(4,1) = 23100929453385.d0/9235860466888.d0
      bD(4,2) = -406327119008035.d0/31390294737638.d0
      bD(4,3) = 58557019738192.d0/4831990238399.d0
      bD(4,4) = -41370260092090.d0/12760377281559.d0
      bD(5,1) = -10008660230968.d0/2813643758547.d0
      bD(5,2) = 136988604098439.d0/4786347071965.d0
      bD(5,3) = -739834933727873.d0/20667726835317.d0
      bD(5,4) = 134755361191288.d0/10036583072473.d0
      bD(6,1) = 4549150143932.d0/17102519410635.d0
      bD(6,2) = -12746474459162.d0/17847296939681.d0
      bD(6,3) = 5543576992611.d0/7464798720697.d0
      bD(6,4) = -94251174957.d0/5756881231433.d0
    
      elseif(icase.eq.6) then     !    The scheme you're using is: ESDIRK436L2SA_C3
    
      ns = 6 

c     The scheme coefficients, A, b, bh, and c  are given by
    
      ai(2,1) = 1.d0/4.d0
      ai(2,2) = 1.d0/4.d0
      ai(3,1) = 1.d0/8.d0
      ai(3,2) = -685497852586.d0/11816340736199.d0
      ai(3,3) = 1.d0/4.d0
      ai(4,1) = 1661069150023.d0/17910800065638.d0
      ai(4,2) = -6217669272415.d0/21853084026797.d0
      ai(4,3) = 2398491653538.d0/4654100100763.d0
      ai(4,4) = 1.d0/4.d0
      ai(5,1) = -291533403691.d0/3304616562838.d0
      ai(5,2) = -12871587167682.d0/13399274912851.d0
      ai(5,3) = 10247774747950.d0/8599160607423.d0
      ai(5,4) = 6078955942878.d0/10469686453393.d0
      ai(5,5) = 1.d0/4.d0
      ai(6,1) = 1189740399107.d0/8132927326164.d0
      ai(6,3) = 1428229024078.d0/7147164224233.d0
      ai(6,4) = 304309715226.d0/589294987153.d0
      ai(6,5) = -814663191941.d0/7240506707015.d0
      ai(6,6) = 1.d0/4.d0
    
      be(1) = 1189740399107.d0/8132927326164.d0
      be(2) = 0.d0/1.d0
      be(3) = 1428229024078.d0/7147164224233.d0
      be(4) = 304309715226.d0/589294987153.d0
      be(5) = -814663191941.d0/7240506707015.d0
      be(6) = 1.d0/4.d0
    
      beh(1) = 355483759106.d0/3424010275803.d0
      beh(2) = -1625665259487.d0/9316172364508.d0
      beh(3) = 1217971653666.d0/3203566609955.d0
      beh(4) = 3544956981031.d0/6728920745444.d0
      beh(5) = -751304069141.d0/22093822668905.d0
      beh(6) = 3037153856030.d0/15364990013649.d0
    
      bi(1) = 1189740399107.d0/8132927326164.d0
      bi(2) = 0.d0/1.d0
      bi(3) = 1428229024078.d0/7147164224233.d0
      bi(4) = 304309715226.d0/589294987153.d0
      bi(5) = -814663191941.d0/7240506707015.d0
      bi(6) = 1.d0/4.d0
    
      bih(1) = 355483759106.d0/3424010275803.d0
      bih(2) = -1625665259487.d0/9316172364508.d0
      bih(3) = 1217971653666.d0/3203566609955.d0
      bih(4) = 3544956981031.d0/6728920745444.d0
      bih(5) = -751304069141.d0/22093822668905.d0
      bih(6) = 3037153856030.d0/15364990013649.d0
    
c   The primary predictor is now given in two parts; eta and B
    
c   U^{(n+1,i+1)} = \sum_{j=0}^{order}eta_{ij}*r^{j} +  
c                 + \sum_{j=1}^{S-1}\sum_{k=0}^{order}B_{ijk}*r^{k}U^{(n,j+1)} 
    
      svpB(2,6,0) = 1.d0/1.d0
      svpB(3,6,0) = 1.d0/1.d0
      svpB(4,6,0) = 1.d0/1.d0
      svpB(5,6,0) = 1.d0/1.d0
      svpB(6,6,0) = 1.d0/1.d0
      svpB(2,2,1) = 622820019924738.d0/5168760291587.d0
      svpB(3,2,1) = 554310012178316.d0/9241757992777.d0
      svpB(4,2,1) = 511637792870717.d0/11860493926638.d0
      svpB(5,2,1) = -452031555857547.d0/9587130357494.d0
      svpB(6,2,1) = 564570624138349.d0/8344116840377.d0
      svpB(2,3,1) = -1232968300451843.d0/7515672190266.d0
      svpB(3,3,1) = -757809040366714.d0/9344449903879.d0
      svpB(4,3,1) = -491053786985383.d0/8851911278036.d0
      svpB(5,3,1) = 114743903081056.d0/1553888606361.d0
      svpB(6,3,1) = -804746700661816.d0/9337785405907.d0
      svpB(2,4,1) = 115040036572611.d0/7322742831637.d0
      svpB(3,4,1) = 150900687268297.d0/21491433528318.d0
      svpB(4,4,1) = 10351943808510.d0/10365195366451.d0
      svpB(5,4,1) = -125432146297267.d0/6318067099543.d0
      svpB(6,4,1) = 233245753037.d0/598263318076.d0
      svpB(2,5,1) = 380181685392696.d0/9323326641643.d0
      svpB(3,5,1) = 276478853415593.d0/13480638955837.d0
      svpB(4,5,1) = 238138283937667.d0/15047609063705.d0
      svpB(5,5,1) = -237136823667953.d0/19249070464867.d0
      svpB(6,5,1) = 273362623341931.d0/10876049968873.d0
      svpB(2,6,1) = -1868655275268133.d0/33100731971985.d0
      svpB(3,6,1) = -428207662509307.d0/15315640595841.d0
      svpB(4,6,1) = -797930134293895.d0/41151745823186.d0
      svpB(5,6,1) = 223203848115757.d0/9102510542238.d0
      svpB(6,6,1) = -265575518177656.d0/8792735266971.d0
      svpB(2,2,2) = -437336324643850.d0/9594031586743.d0
      svpB(3,2,2) = -70097533931319.d0/3825996964615.d0
      svpB(4,2,2) = -1267202575112647.d0/21125094304614.d0
      svpB(5,2,2) = -1225642649490109.d0/7092735396475.d0
      svpB(6,2,2) = -2118769417914809.d0/11620084803757.d0
      svpB(2,3,2) = 411744113858521.d0/5826318059332.d0
      svpB(3,3,2) = 365272875006371.d0/12859986157972.d0
      svpB(4,3,2) = 740105976733313.d0/7958435806330.d0
      svpB(5,3,2) = 1710865422940590.d0/6386265048451.d0
      svpB(6,3,2) = 1821712951478755.d0/6444463145238.d0
      svpB(2,4,2) = -130849489404241.d0/8293727189849.d0
      svpB(3,4,2) = -8654273268544.d0/1364786149387.d0
      svpB(4,4,2) = -242899550538850.d0/11699611510863.d0
      svpB(5,4,2) = -1686960590363127.d0/28206353899783.d0
      svpB(6,4,2) = -275645062606351.d0/4367852256999.d0
      svpB(2,5,2) = -26014490620266.d0/2270389192691.d0
      svpB(3,5,2) = -92380161945919.d0/20059494841962.d0
      svpB(4,5,2) = -78165962976284.d0/5184055988269.d0
      svpB(5,5,2) = -911447759642040.d0/20983651565183.d0
      svpB(6,5,2) = -271623058612420.d0/5926409109457.d0
      svpB(2,6,2) = 296212380325894.d0/14383116823771.d0
      svpB(3,6,2) = 99652528136983.d0/12039110797371.d0
      svpB(4,6,2) = 163744033798463.d0/6042022012957.d0
      svpB(5,6,2) = 859699696407263.d0/11011864919207.d0
      svpB(6,6,2) = 519373619423915.d0/6304776521444.d0
      svpB(2,2,3) = -354376503376919.d0/4755622841317.d0
      svpB(3,2,3) = -148586800720349.d0/11738116330399.d0
      svpB(4,2,3) = -1111158556971113.d0/14817000804303.d0
      svpB(5,2,3) = -11619872834237246.d0/31690680876151.d0
      svpB(6,2,3) = -5344388860004255.d0/13447507890614.d0
      svpB(2,3,3) = 629356194069592.d0/5826220487643.d0
      svpB(3,3,3) = 298269310476987.d0/16254539388275.d0
      svpB(4,3,3) = 778828600705872.d0/7164305232749.d0
      svpB(5,3,3) = 1998945125618970.d0/3760789162081.d0
      svpB(6,3,3) = 10559316659466005.d0/18328523945363.d0
      svpB(2,4,3) = -505282494524824.d0/29276466003265.d0
      svpB(3,4,3) = -53849506663771.d0/18367151981511.d0
      svpB(4,4,3) = -216916877892482.d0/12488765978373.d0
      svpB(5,4,3) = -245619328938945.d0/2892242134516.d0
      svpB(6,4,3) = -1129775333394881.d0/12273763747706.d0
      svpB(2,5,3) = -56726587475243.d0/2726620188956.d0
      svpB(3,5,3) = -40640686123843.d0/11499387555192.d0
      svpB(4,5,3) = -202073832250535.d0/9651389644059.d0
      svpB(5,5,3) = -521501787670531.d0/5094264596326.d0
      svpB(6,5,3) = -867110693831401.d0/7814726840695.d0
      svpB(2,6,3) = 222483535427685.d0/6707367209203.d0
      svpB(3,6,3) = 90842573931389.d0/16122015300718.d0
      svpB(4,6,3) = 165593907793787.d0/4960670087245.d0
      svpB(5,6,3) = 562554469909061.d0/3446723897203.d0
      svpB(6,6,3) = 28162148032093445.d0/159192005956814.d0
      svpB(2,2,4) = -333621467255574.d0/9519489212837.d0
      svpB(3,2,4) = 57784669570116.d0/50097917169145.d0
      svpB(4,2,4) = -45655294725773.d0/2881038481565.d0
      svpB(5,2,4) = -2018957493431648.d0/10259728486365.d0
      svpB(6,2,4) = -3241586627327091.d0/14794036811692.d0
      svpB(2,3,4) = 233494073513827.d0/4755461095051.d0
      svpB(3,3,4) = -15828999653536.d0/9858339408057.d0
      svpB(4,3,4) = 272252757413006.d0/12217322603281.d0
      svpB(5,3,4) = 1901038271058864.d0/6883102823063.d0
      svpB(6,3,4) = 2076016900500156.d0/6755512592609.d0
      svpB(2,4,4) = -39658727119123.d0/6060606764246.d0
      svpB(3,4,4) = 330991980931.d0/1652077214318.d0
      svpB(4,4,4) = -21254725291381.d0/6903634083395.d0
      svpB(5,4,4) = -468190420181490.d0/12499724118239.d0
      svpB(6,4,4) = -504567265698481.d0/12192419691664.d0
      svpB(2,5,4) = -83569348320103.d0/8029306027059.d0
      svpB(3,5,4) = 3317317165888.d0/8330836751893.d0
      svpB(4,5,4) = -123017992500893.d0/28865208704375.d0
      svpB(5,5,4) = -521323946017594.d0/9342986474541.d0
      svpB(6,5,4) = -630004088789632.d0/9948639444627.d0
      svpB(2,6,4) = 173088962927553.d0/10924247728592.d0
      svpB(3,6,4) = -8957531605202.d0/15706452772511.d0
      svpB(4,6,4) = 18444995025624.d0/2722751553739.d0
      svpB(5,6,4) = 317602167003383.d0/3665447774886.d0
      svpB(6,6,4) = 1448125899420860.d0/14848157267837.d0
    
c   U^{(n+1,i+1)} = \sum_{j=0}^{order}eta_{ij}*r^{j} +  
c                 + \sum_{j=1}^{S-1}\sum_{k=0}^{order}B_{ijk}*r^{k}U^{(n,j+1)} 
    
      svpB(2,1,1) = 461115308777364.d0/10594959869233.d0
      svpB(3,1,1) = 120863241740264.d0/5609422122881.d0
      svpB(4,1,1) = 14881971164725.d0/998662116427.d0
      svpB(5,1,1) = -736513858262371.d0/38678288321618.d0
      svpB(6,1,1) = 178548621664378.d0/7695874125463.d0
      svpB(2,1,2) = -203166855506287.d0/11014827806010.d0
      svpB(3,1,2) = -90075104801316.d0/12150269694955.d0
      svpB(4,1,2) = -59902840418731.d0/2467966643890.d0
      svpB(5,1,2) = -210591756941811.d0/3011830460986.d0
      svpB(6,1,2) = -994996663045733.d0/13486103434160.d0
      svpB(2,1,3) = -282503547445894.d0/9874197766131.d0
      svpB(3,1,3) = -65517653660367.d0/13480684457423.d0
      svpB(4,1,3) = -250655732867129.d0/8705572592747.d0
      svpB(5,1,3) = -1173627438456163.d0/8336732743085.d0
      svpB(6,1,3) = -613927071775499.d0/4023425574320.d0
      svpB(2,1,4) = -171382260755507.d0/13237434601144.d0
      svpB(3,1,4) = 30311478566405.d0/71493762901244.d0
      svpB(4,1,4) = -52285338483733.d0/8905463046620.d0
      svpB(5,1,4) = -2039111715871749.d0/28010721890125.d0
      svpB(6,1,4) = -5002482991369131.d0/61749643354804.d0
    
c  The secondary predictor is used well into the next step
c  It uses only information from the new step 
    
c   Z_{0}^{(n+1,i)} = \sum_{j=2}^{i-1} \alpha_{ij} Z_{0}^{(n+1,j)} 
    
      alpha(3,2) = 7141075053842.d0/11263976658481.d0
      alpha(4,2) = 14534263011338.d0/9037123749793.d0
      alpha(4,3) = -12145565496278.d0/16697561406817.d0
      alpha(5,2) = 47814865418801.d0/10432170648376.d0
      alpha(5,3) = -53336429163275.d0/7775119750182.d0
      alpha(5,4) = 47984328140430.d0/32141113519727.d0
      alpha(6,2) = 14097636665547.d0/5702313947045.d0
      alpha(6,3) = -73974254437330.d0/21204619481989.d0
      alpha(6,4) = 4249616235173.d0/8536665442503.d0
      alpha(6,5) = 2388038697443.d0/3979545713606.d0
    
c   This secondary predictor uses information from the new and old step 
    
    
c   U^{(n+1,3)} = al3_{1}*U^{(n  ,3)} + al3_{2}*U^{(n  ,4)} +
c               + al3_{3}*U^{(n  ,5)} + al3_{4}*U^{(n  ,6)} +
c               + al3_{5}*U^{(n+1,2)}
    
    
c   al3_{i} = \sum_{j=1}^{2*(order)} al3N(i,j)*r^{j} / 
c             \sum_{j=1}^{2*(order)} al3D(i,j)*r^{j} 
    
    
      al3N(1,1) = -21372682853.d0/12646589952648.d0
      al3N(2,1) = 73219228104.d0/10287737095451.d0
      al3N(3,1) = -10975443751393.d0/15233235636567.d0
      al3N(4,1) = 7731402997495.d0/5887711301337.d0
      al3N(5,1) = 1657092233154.d0/4122901604639.d0
      al3N(1,2) = 1467118481979.d0/10292593534099.d0
      al3N(2,2) = -47242008094463.d0/80611965953417.d0
      al3N(3,2) = 12929215662019.d0/9379534655471.d0
      al3N(4,2) = -4832431147767.d0/11361861064642.d0
      al3N(5,2) = 1575038286103.d0/15102813431212.d0
      al3N(1,3) = 3717603395893.d0/6735962588790.d0
      al3N(2,3) = -14495196667317.d0/10143716478025.d0
      al3N(3,3) = 15891644425219.d0/4681367554254.d0
      al3N(4,3) = -17649892325588.d0/6871019832861.d0
      al3N(5,3) = 27617951645.d0/9221106656096.d0
      al3N(1,4) = 3120008636498.d0/13021228562315.d0
      al3N(2,4) = -10881333329982.d0/25677079971559.d0
      al3N(3,4) = 36492656803013.d0/28588545795884.d0
      al3N(4,4) = -9405543489288.d0/8610678683155.d0
      al3N(1,5) = 161446249924.d0/7337993530617.d0
      al3N(2,5) = -1577996795474.d0/17413507893585.d0
      al3N(3,5) = 2620300858066.d0/2939888963453.d0
      al3N(4,5) = -35251777140231.d0/42850199399204.d0
      al3N(1,6) = 170647587001.d0/9367963012150.d0
      al3N(2,6) = -89108781280.d0/1860488356307.d0
      al3N(3,6) = 1578616204324.d0/5239983264561.d0
      al3N(4,6) = -3608954755977.d0/13288524621397.d0
    
    
      al3D(1,1) = 1.d0/1.d0
      al3D(2,1) = 1.d0/1.d0
      al3D(3,1) = 1.d0/1.d0
      al3D(4,1) = 1.d0/1.d0
      al3D(5,1) = 1.d0/1.d0
      al3D(1,2) = 4129438051769.d0/6726389405720.d0
      al3D(2,2) = 4129438051769.d0/6726389405720.d0
      al3D(3,2) = 4129438051769.d0/6726389405720.d0
      al3D(4,2) = 4129438051769.d0/6726389405720.d0
      al3D(5,2) = 4129438051769.d0/6726389405720.d0
      al3D(1,3) = -472359906729.d0/9806186422679.d0
      al3D(2,3) = -472359906729.d0/9806186422679.d0
      al3D(3,3) = -472359906729.d0/9806186422679.d0
      al3D(4,3) = -472359906729.d0/9806186422679.d0
      al3D(5,3) = -472359906729.d0/9806186422679.d0
    
    
c   This secondary predictor uses information from the new and old step 
    
    
c   U^{(n+1,4)} = al4_{1}*U^{(n  ,4)} + al4_{2}*U^{(n  ,5)} +
c               + al4_{3}*U^{(n  ,6)} + al4_{4}*U^{(n+1,2)} +
c               + al4_{5}*U^{(n+1,3)}
    
    
c   al4_{i} = \sum_{j=1}^{2*(order)} al4N(i,j)*r^{j} / 
c             \sum_{j=1}^{2*(order)} al4D(i,j)*r^{j} 
    
    
      al4N(2,1) = -1.d0/476371321430964598863187.d0
      al4N(3,1) = 2774235469737.d0/2652524679092.d0
      al4N(4,1) = 15320879859895.d0/6866936769868.d0
      al4N(5,1) = -40311964614133.d0/17704033663015.d0
      al4N(1,2) = 45815432410.d0/28038048561109.d0
      al4N(2,2) = -61923654343518.d0/9090295058315.d0
      al4N(3,2) = 58858061204038.d0/7740133013645.d0
      al4N(4,2) = 17422177602561.d0/4286895422063.d0
      al4N(5,2) = -4079631132238.d0/4137529424393.d0
      al4N(1,3) = 739239120159.d0/9491306572688.d0
      al4N(2,3) = -263587989251166.d0/13065427253675.d0
      al4N(3,3) = 228287003142797.d0/12522690611188.d0
      al4N(4,3) = 14154435121833.d0/10844357951093.d0
      al4N(5,3) = 55856389974937.d0/24908555749640.d0
      al4N(1,4) = 45406042302539.d0/22363479551451.d0
      al4N(2,4) = -193118717137693.d0/10968014258793.d0
      al4N(3,4) = 144267812217394.d0/9423097669611.d0
      al4N(4,4) = 231275906774.d0/4099573397957.d0
      al4N(5,4) = 2157354458899.d0/5910609406200.d0
      al4N(1,5) = -6340502891.d0/12598459156872.d0
      al4N(2,5) = 15468318523054.d0/7372626661679.d0
      al4N(3,5) = -20071085112568.d0/8163649250139.d0
      al4N(4,5) = 724093733683.d0/13189683331319.d0
      al4N(5,5) = 5013706131028.d0/11554501468495.d0
      al4N(1,6) = -148148999574.d0/10399376825737.d0
      al4N(2,6) = 3769085408134.d0/1021413303325.d0
      al4N(3,6) = -88713160656741.d0/24134231488411.d0
      al4N(1,7) = -715600761283.d0/7985585513433.d0
      al4N(2,7) = 3698352211610.d0/2564322889689.d0
      al4N(3,7) = -14492752736464.d0/10714563306431.d0
    
    
      al4D(1,1) = 1.d0/1.d0
      al4D(2,1) = 1.d0/1.d0
      al4D(3,1) = 1.d0/1.d0
      al4D(4,1) = 1.d0/1.d0
      al4D(5,1) = 1.d0/1.d0
      al4D(1,2) = 53961028886899.d0/13936611634887.d0
      al4D(2,2) = 53961028886899.d0/13936611634887.d0
      al4D(3,2) = 53961028886899.d0/13936611634887.d0
      al4D(4,2) = 53961028886899.d0/13936611634887.d0
      al4D(5,2) = 53961028886899.d0/13936611634887.d0
      al4D(1,3) = 15033076749071.d0/8943024210186.d0
      al4D(2,3) = 15033076749071.d0/8943024210186.d0
      al4D(3,3) = 15033076749071.d0/8943024210186.d0
      al4D(4,3) = 15033076749071.d0/8943024210186.d0
      al4D(5,3) = 15033076749071.d0/8943024210186.d0
      al4D(1,4) = 496642385686.d0/3217610793455.d0
      al4D(2,4) = 496642385686.d0/3217610793455.d0
      al4D(3,4) = 496642385686.d0/3217610793455.d0
      al4D(4,4) = 496642385686.d0/3217610793455.d0
      al4D(5,4) = 496642385686.d0/3217610793455.d0
      al4D(1,5) = 671999293417.d0/5258407006027.d0
      al4D(2,5) = 671999293417.d0/5258407006027.d0
      al4D(3,5) = 671999293417.d0/5258407006027.d0
      al4D(4,5) = 671999293417.d0/5258407006027.d0
      al4D(5,5) = 671999293417.d0/5258407006027.d0
    
c   Dense output 
    
c   bd_{i} = \sum_{j=1}^{dense_order}  bdense_{ij}*x^{j} 
    
c   U[t^{n} + x*(Delta t)] = U^{(n)}  
c           + (Delta t)\sum_{i=1}^{s} bd_{i}*F^{(i)} 
    
    
      bD(1,1) = 6947827441573.d0/8903452997917.d0
      bD(1,2) = -15151522803474.d0/10298541125287.d0
      bD(1,3) = 11238834868769.d0/9471848008518.d0
      bD(1,4) = -885193820293.d0/2533565647497.d0
      bD(2,1) = -1.d0/801917298145728384576274.d0
      bD(2,2) = 1.d0/120632267811045841319720.d0
      bD(2,3) = -1.d0/79632010182174796171585.d0
      bD(2,4) = 1.d0/180231720199587795006835.d0
      bD(3,1) = 11186832770183.d0/10494404909745.d0
      bD(3,2) = -5266192794660.d0/2620338527339.d0
      bD(3,3) = 5391770640018.d0/3326488786783.d0
      bD(3,4) = -3385207121193.d0/7092841553977.d0
      bD(4,1) = -10936654328996.d0/7060593836741.d0
      bD(4,2) = 341761990967577.d0/43280765257963.d0
      bD(4,3) = -74148770233127.d0/8165902326448.d0
      bD(4,4) = 15376862021002.d0/4732414065377.d0
      bD(5,1) = 44266881533842.d0/7711412251119.d0
      bD(5,2) = -331297217760877.d0/8577657277705.d0
      bD(5,3) = 545047166648793.d0/9148896881783.d0
      bD(5,4) = -262220737224484.d0/9782583379485.d0
      bD(6,1) = -47787731276353.d0/9485834337250.d0
      bD(6,2) = 2223837670594285.d0/65009574471218.d0
      bD(6,3) = -429336256813307.d0/8054741903376.d0
      bD(6,4) = 132807060383020.d0/5446874804559.d0
    
       endif                                                  !  end icase loop
c
       do i=1,ns
         bi(i)  = be(i)
         bih(i) = beh(i)
       enddo

       do i = 1,ns
         do j = 1,ns
           cE(i) = cE(i) + aE(i,j)
           cI(i) = cI(i) + aI(i,j)
         enddo
       enddo
 
       eps = 1.0e-7
       if(beh(1).lt.eps.and.bih(1).lt.eps)then
       do i = 1,ns
         beh(i) = 1./ns
         bih(i) = 1./ns
       enddo
       endif

       return
       end

      SUBROUTINE fit(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q)
      INTEGER mwt ,ndata
      REAL a,b,chi2,q,siga,sigb,sig(ndata),x(ndata),y(ndata)
      INTEGER i
      REAL sigdat,ss,st2,sx,sxoss,sy,t,wt,gammq
      sx=0.
      sy=0.
      st2=0. 
      b=0.
      if(mwt.ne.0) then 
        ss=0.
        do i=1,ndata
          wt=1./(sig(i)**2)
          ss=ss+wt
          sx=sx+x(i)*wt
          sy=sy+y(i)*wt
        enddo
      else
        do i=1,ndata
          sx=sx+x(i)
          sy=sy+y(i)
        enddo 
        ss=float(ndata)
      endif
      sxoss=sx/ss
      if(mwt.ne.0) then
        do i=1,ndata
          t=(x(i)-sxoss)/sig(i)
          st2=st2+t*t
          b=b+t*y(i)/sig(i)
        enddo
      else
        do i=1,ndata
          t=x(i)-sxoss
          st2=st2+t*t
          b=b+t*y(i)
        enddo 
      endif
      b=b/st2
      a=(sy-sx*b)/ss
      siga=sqrt((1.+sx*sx/(ss*st2))/ss)
      sigb=sqrt(1./st2)
      chi2=0.
      if(mwt.eq.0) then
        do i=1,ndata
          chi2=chi2+(y(i)-a-b*x(i))**2
        enddo 
        q=1.
        sigdat=sqrt(chi2/(ndata-2) )
        siga=siga*sigdat
        sigb=sigb*sigdat
      else
        do  i=1,ndata
          chi2=chi2+((y(i)-a-b*x(i))/sig(i))**2
        enddo 
        q=gammq(0.5*(ndata-2),0.5*chi2)
      endif
      return
      end

      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER itkax
      REAL a,gamser,gln,x,eps
      PARAMETER (itkax=100,eps=3.e-7)
c     USES gammln
c     Returns the incomplete gamma function P(a, x) evaluated by its series 
c     representation as gamser. Also returns In r( a) as gln.
      INTEGER n
      REAL ap,del,sum,gammln
      gln=gammln (A)
      if(x.le.0.)then
      if(x.lt.0.)pause 'x < 0 in gser'
      gamser=0.
      return
      endif
      ap=a
      sum=1./a
      del=sum
      do n=1,itkax
      ap=ap+1.
      del=del*x/ap
      sum=sum+del
      if(abs(del).lt.abs(sum)*eps)goto 1
      enddo 
      pause 'a too large, itkax too small in gser'
    1 gamser=sum*exp(-x+a*log(x)-gln)
      return
      end

      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER itkax
      REAL a,gammcf,gln,x,eps,FPMIN
      PARAMETER (itkax=100,eps=3.e-7,FPMIN=1.e-30)
c     USES gammln
      INTEGER i
      REAL an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1. -a 
      c=1. /FPMIN 
      d=1./b
      h=d
      do i=1, itkax 
      an=-i*(i-a)
      b=b+2.
      d=an*d+b
      if(abs(d).lt.FPMIN)d=FPMIN
      c=b+an/c
      if(abs(c).lt.FPMIN)c=FPMIN
      d=1./d
      del=d*c
      h=h*del
      if(abs(del-1.).lt.EPS)goto 1
      enddo
      pause 'a too large, ITMAX too small in gcf'
    1 gammcf=exp(-x+a*log(x)-gln)*h
      return
      end

      FUNCTION gammq(a,x)
      REAL a,gammq,x
c     USES gcf,gser
c     Returns the incomplete gamma function Q(a, x) = 1 -P(a, x).
      REAL gammcf ,gamser ,gln
      if(x.lt.0. .or.a.le.0.)pause 'bad arguments in gammq'
      if (x .lt .a+1. ) then 
      call gser(gamser,a,x,gln)
      gammq=1. -gamser 
      call gcf(gammcf,a,x,gln)
      endif
      return
      end

      FUNCTION gammln(xx)
      REAL gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     &24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     &-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
      y=y+1.d0
      ser=ser+cof(j)/y
      enddo 
      gammln=tmp+log(stp*ser/x)
      return
      end
