c  program obliczajacy szum w przedziale bmin - bmax w (c/d) dla FT 
c   *.trf

c uzycie:
c
c s2n_cd  fr1  fr2

      double precision bmin, bmax, t, flux, flux1,aver,a, sigma4
      character*80 inname, nf1,nf2
 
      call getarg(1,inname)
      call getarg(2,nf1)
      call getarg(3,nf2)
      
      read(nf1,*) bmin
      read(nf2,*) bmax      
      
c      bmin=0.004d+0   ! 0.864 c/d to 10.0 uHz
c      bmax=0.02d+0       ! 6.0   c/d to 69.5 uHz
      i=0
      
      a=0.d+0
      flux1=0.d+0
c     otworz i wczytaj kol. 1, 2  z FT (*.tdat.trf) 
      open(1,file=inname,err=135,status='old')
      do j=1, 5000000
        if(j.eq.5000000) then
          write(*,*) ' declared MATRIX too small!! EXIT'
          goto 300
        endif
        read(1,*,end=53) t, flux
        if(t.ge.bmin .and. t.le.bmax) then
            flux1=flux1+flux
            i=i+1
        endif 
      enddo
53    close(1)

      a=i-1.d+0
      aver=flux1/a
      sigma4=aver*4.d+0

c      write(*,*) 'aver ',aver
c      write(*,*) 'sigma4 ', sigma4
       
c      write(*,*)
      write(*,'("S/N between: ",f10.5," - ",f10.5," c/d")') bmin,bmax
c      write(*,*)
      write(*,'("average: ",f16.9,1x," 4sig: ",f16.9," 5sig: ",f16.9)') 
     # aver, sigma4, sigma4+aver

c      write(*,*)
c      write(*,'(" average in sec: ",f16.9,1x," 4sigma ",f16.9)')
c     #  aver*86400d+0, sigma4*86400.d+0

      write(*,*)
      goto 300
      
135   write(*,*) ' error opening file: ', inname
 
300   end
