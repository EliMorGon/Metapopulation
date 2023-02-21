module comun
implicit none
integer solopt,nind,ntotal,npob,estruc,ncrom,nm,nh
double precision particion
allocatable solopt(:,:,:),ntotal(:,:),estruc(:,:)
dimension particion(3,2)
end module comun

module decision
implicit none
integer nmigra
double precision fped,lambda
double precision fgeno
allocatable fgeno(:,:),fped(:,:)
end module decision


program meta_gen

! Gestión automática de metapoblaciones                          *
! Seleccion por minimo parentesco genealogico o genomico         *
! Numero de machos y hembras iguales o diferentes                *
! Numero de marcadores genotipados variable                      *
! Contribuciones diferenciales y apareamiento al azar            *

use comun
use decision

implicit none
integer*1  pop,marcadores,unobserved
integer popt1,o,geneal,pos1,pos2,ran1,n_SNP_total,&
        matriz,ngen,nrep,i,j,t,k,l,m,movil,&
        pos,alelo,cross,sobant,pad,n,igs,ngs,&
        irep,igen,nsob,aux,total,jcr,&
        ipedcam,isex,marct1,totalmar,&
        migrado,kkkk,contB,contN,contN1,contB1,kk,totalneut
integer genom_length,n_neut_total,n_neut,n_SNP,n_pos,genom_map,SNP_pos,neut_pos,&
        genom_map_neut,tfreq
double precision long,u,results,inter,sk,sk2,param1,param2,interval,&
     fem,conm,homgenom,ibdmar,homar,divgenmar,pOBSS,pOBSSpop,pOBSSpop2,long0,&
     dival,divgen,fpedt1,fpm,conped,divalmar,fmol,pOBSSpop1,pOBSNpop,&
     pOBSNpop1,divgenobs,divalobs,time1,time2,time3,fmar,conmar
double precision divgenobsn,divalobsn,homgenomn,ibdmarn,homarn,divgenneut,&
        divalneut,fmoln,fmarn,conmarn,results2,results3,&
        homgenoms,ibdmars,homars,homgenomns,ibdmarns,homarns
character arep*2
character(len=11) :: a
integer*8 total_genom_length
allocatable pop(:,:,:,:),popt1(:,:,:,:),geneal(:,:),movil(:,:,:),&
            results(:,:),marct1(:,:,:,:),totalmar(:,:),&
            cross(:),aux(:),fmol(:,:),param1(:),param2(:),fmar(:,:),&
            fpedt1(:,:),marcadores(:,:,:,:),fem(:,:),fpm(:,:),&
            ipedcam(:,:),total(:,:),conmar(:),pOBSSpop1(:,:),conped(:),conm(:),&
            contB(:),contB1(:),pOBSS(:),pOBSSpop(:,:),pOBSNpop(:,:),pOBSSpop2(:,:,:),&
            interval(:),contN(:),pOBSNpop1(:,:),contN1(:)
allocatable unobserved(:,:,:,:),results2(:,:),results3(:,:,:),&
            conmarn(:),fmarn(:,:),fmoln(:,:),totalneut(:,:)
allocatable genom_length(:),n_neut(:),n_SNP(:),n_pos(:),genom_map(:,:),&
            SNP_pos(:,:),neut_pos(:,:),genom_map_neut(:,:),homgenoms(:),&
            ibdmars(:),homars(:),homgenomns(:),ibdmarns(:),homarns(:)
dimension sk(2,2),sk2(2,2),migrado(2)

! Lee especificaciones de la simulacion *
open(5,FILE="datos_metapopulation.txt")
! N. generaciones, N. generaciones sin control, N. repeticiones *
read(5,*) ngen,ngs,nrep
! N. poblaciones *
read(5,*) npob
allocate(ntotal(npob,2),estruc(npob,2))
! N. machos y hembras por población *
do i=1,npob
  read(5,*) (ntotal(i,j),j=1,2)
enddo
nind=sum(ntotal)
nm=sum(ntotal(:,1))
nh=sum(ntotal(:,2))
estruc(1,1)=1
estruc(1,2)=estruc(1,1)+ntotal(1,1)
do i=2,npob
  estruc(i,1)=sum(ntotal(1:i-1,:))+1
  estruc(i,2)=estruc(i,1)+ntotal(i,1)
enddo
! N. cromosomas! All chromosomes of the same length !
read(5,*) ncrom,total_genom_length
  	print*,total_genom_length


! Peso de la consanguinidad y numero maximo de migrantes *
read(5,*) lambda,nmigra
nmigra=nmigra*2

! matriz de parentesco a usar *
read(5,*) matriz

!frecuencias que usar
read(5,*) tfreq


close(5)

         open(UNIT=84,FILE='matriz.txt')

! Reserva espacio para matrices diversas *


allocate (results(0:ngen,41),cross(1),aux(nind),ipedcam(nind,2))
allocate (results2(0:ngen,41),results3(npob,0:ngen,12))
allocate (fped(nind,nind),fpedt1(nind,nind))
allocate (fmol(nind,nind),fgeno(nind,nind))
allocate (fmoln(nind,nind))
allocate (solopt(nind,npob,2),movil(npob,npob,ngen))
allocate (fem(npob+1,npob+1),fpm(npob+1,npob+1),fmar(npob+1,npob+1))
allocate (fmarn(npob+1,npob+1))
allocate (conm(npob+1),conped(npob+1),conmar(npob+1))
allocate (geneal(nind,2))
allocate(homgenoms(npob),ibdmars(npob),homars(npob),homgenomns(npob))
allocate(ibdmarns(npob),homarns(npob),conmarn(npob+1))


results2=0
results3=0
results=0
movil=0


call cpu_time(time1)
time2=time1

allocate(n_neut(ncrom),n_SNP(ncrom),n_pos(ncrom),genom_length(ncrom),interval(ncrom))  

! Bucle para repeticiones *
do irep=1,nrep

    open(UNIT=35,FILE='repl.txt')
    write(35,*) irep
    close(35)
    
        write(arep,'(i2)') irep
    arep=adjustl(arep)

      print*,'Replica',irep
      call setup(10*(irep+0))
      
      call SYSTEM('./run.cor')
      
      


n_neut=0
n_SNP=0
n_pos=0
genom_length=0
interval=0
   
        
! Read  a file the genome maps characteristics !


open(UNIT=35,FILE='newmaps.txt',status='old')
read(35,*)
do i=1,ncrom
 read(35,*) a,jcr,genom_length(i),n_SNP(i),n_neut(i)
enddo
allocate(SNP_pos(ncrom,maxval(n_SNP)))
allocate(genom_map(ncrom,maxval(n_SNP)))
SNP_pos=0
genom_map=0


n_pos=n_neut+n_SNP
print*,n_pos



read(35,*)
do i=1,ncrom
  read(35,*)
  do j=1,n_SNP(i)
    read(35,*) SNP_pos(i,j),genom_map(i,j)
  enddo
enddo

close(35)
!close(352)

n_SNP_total=sum(n_SNP)
n_neut_total=sum(n_neut)

open(UNIT=35,FILE='newmapsneut.txt',status='old')


allocate(neut_pos(ncrom,maxval(n_neut)))
allocate(genom_map_neut(ncrom,maxval(n_neut)))

neut_pos=0
genom_map_neut=0



do i=1,ncrom 
  read(35,*)
  do j=1,n_neut(i)
    read(35,*) neut_pos(i,j),genom_map_neut(i,j)
  enddo
enddo

close(35)


allocate (totalmar(ncrom,maxval(n_SNP)),total(ncrom,maxval(n_pos)))
allocate (totalneut(ncrom,maxval(n_neut)))


do j=1,ncrom
    k=0
    do i=1,n_pos(j)
        k=k+1
        total(j,i)=k
    enddo
enddo
do j=1,ncrom
    k=0    
    do i=1,n_SNP(j)
        k=k+1
        totalmar(j,i)=k
    enddo
enddo

do j=1,ncrom
    k=0    
    do i=1,n_neut(j)
        k=k+1
        totalneut(j,i)=k
    enddo
enddo


allocate (pop(2,maxval(n_pos),ncrom,nind),popt1(2,maxval(n_pos),ncrom,nind))
allocate (marcadores(2,maxval(n_SNP),ncrom,nind))
allocate (unobserved(2,maxval(n_neut),ncrom,nind))
allocate (marct1(2,maxval(n_SNP),ncrom,nind))

! Genotypes of the base population !
! Read the initial genotypes of the base population !

open(UNIT=25,FILE='base.txt',status='old')
    do i=1,nind
        !print*,i
        do j=1,ncrom
            read(25,*) k,l,((pop(k,l,j,i),k=1,2),l=1,n_pos(j))
        enddo
    enddo


    close(25)
      
      do i=1,nind
          do j=1,ncrom
              do l=1,n_SNP(j)
                      marcadores(:,l,j,i)=pop(:,SNP_pos(j,l),j,i)  
              enddo
          enddo
      enddo
      
      
      do i=1,nind
          do j=1,ncrom
              do l=1,n_neut(j)
                      unobserved(:,l,j,i)=pop(:,neut_pos(j,l),j,i)  
              enddo
          enddo
      enddo        
          
          

    
    
    ! Inicializa las matrices de parentescos *
      fmol=0
      fped=0
      do i=1,nind
        fped(i,i)=.5
      enddo

    ! Bucle para generaciones sin control*
      do igs=1,ngs

    ! Elige padres al azar y los aparea *
        ran1=sum(ntotal(1,:))
	    do i=1,ran1
          call mother(u)
          geneal(i,1)=min(int(u*ntotal(1,1))+1,ntotal(1,1))
	      call mother(u)
          geneal(i,2)=min(int(u*ntotal(1,2))+ntotal(1,1)+1,ran1)
	    enddo

        do i=1+ran1,nind
          call mother(u)
          geneal(i,1)=min(int(u*(nm-ntotal(1,1)))+ran1+1,nm+ran1-ntotal(1,1))
	      call mother(u)
          geneal(i,2)=min(int(u*(nh-ntotal(1,2)))+nm+ran1-ntotal(1,1)+1,nind)
	    enddo


    ! Genera nueva poblacion *

        print*,'Empieza a generar los hijos',igs
    
    do j=1,ncrom
        long0=dble(genom_length(j))/100000000
        long=long0
        do i=1,nind
            do k=1,2

                 
            interval=dble(genom_length)/dble(n_pos)
    ! Decide el numero de sobrecruzamientos y su posicion *
    ! y los ordena                                        *
            
            
    ! For small chromosomes (less than 1 Morgan) !
            if(long.lt.1) then
	            call mother(u)
            
	            if(u.lt.long) then
	                long=1
	            else
		            long=0
                endif
            else
                long=long0
            
            endif
            
              nsob=0
              do n=1,anint(long)
    
                call mother(u)
                !print*,u,'Entra en Poisson'
                call poisson(u*.9999,l)
                nsob=nsob+l
              enddo
                
              deallocate (cross)
              allocate (cross(nsob))
              cross=0
              do n=1,nsob
                call mother(u)
                cross(n)=min(int(u*genom_length(j)+1),genom_length(j))
                
		        cross(n)=anint(cross(n)/interval(j))
              enddo

              if(nsob.gt.1) call ordena(cross,nsob)

    ! Copia las partes de los homologos correspondientes *
              call mother(u)
              pad=min(int(u*2+1),2)
              sobant=0
              pos=1
              do n=1,nsob
                popt1(k,sobant+1:cross(n),j,i)=pop(pad,sobant+1:cross(n),j,geneal(i,k))
                sobant=cross(n)
                do l=pos,n_SNP(j)
                  if(SNP_pos(j,l).ge.cross(n)) exit
	              marct1(k,l,j,i)=marcadores(pad,l,j,geneal(i,k))
	              pos=pos+1
	            enddo
                select case(pad)
                  case(1)
                    pad=2
                  case(2)
                    pad=1
                end select
              enddo
              popt1(k,sobant+1:n_pos(j),j,i)=pop(pad,sobant+1:n_pos(j),j,geneal(i,k))
	            if(pos.le.n_SNP(j)) then
                   marct1(k,pos:n_SNP(j),j,i)=marcadores(pad,pos:n_SNP(j),j,geneal(i,k))
	            endif
            enddo
          enddo
        enddo
        pop=popt1
        popt1=0
	    marcadores=marct1
	    marct1=0
        unobserved=0

      do i=1,nind
          do j=1,ncrom
              do l=1,n_neut(j)
                      unobserved(:,l,j,i)=pop(:,neut_pos(j,l),j,i)  
              enddo
          enddo
      enddo
      
    ! Calcula parentesco genealogico *
        do i=1,nind
          fpedt1(i,i)=.5*(1+fped(geneal(i,1),geneal(i,2)))
          do j=i+1,nind
            fpedt1(i,j)=.25*(fped(geneal(i,1),geneal(j,1))+&
                             fped(geneal(i,1),geneal(j,2))+&
                             fped(geneal(i,2),geneal(j,1))+&
                             fped(geneal(i,2),geneal(j,2)))
            fpedt1(j,i)=fpedt1(i,j)
          enddo
        enddo
        fped=fpedt1

 

    ! Cierra bucle de generaciones sin control *
      enddo



    ! Calcula valores de partida

    ! Calcula heterocigosidad real *
      call hetero(pop,n_pos,n_SNP,SNP_pos,marcadores,homgenom,ibdmar,homar)
    call hetero(pop,n_pos,n_neut,neut_pos,unobserved,homgenomn,ibdmarn,homarn)
    call hetero2(pop,n_pos,n_SNP,SNP_pos,marcadores,homgenoms,ibdmars,homars)
    call hetero2(pop,n_pos,n_neut,neut_pos,unobserved,homgenomns,ibdmarns,homarns)

    !   print*,'homocigosidad genomica,pos. marcadores, marcadores'
    !	print*,homgenom,ibdmar,homar

    ! Calcula diversidad alelica y genica genomica *
      call diversid(pop,n_pos,n_pos,divgen,dival,total)
      

    !   print*,'diversidades genomicas'
    !	print*,dival,divgen

    ! Calcula diversidad alelica y genica   *
    ! para las posiciones de los marcadores *
      call diversid(pop,n_SNP,n_pos,divgenobs,divalobs,SNP_pos)
    call diversid(pop,n_neut,n_pos,divgenobsn,divalobsn,neut_pos)

    !   print*,'diversidades pos. marcador'
    !	print*,divalobs,divgenobs
                 
    ! Calcula diversidad alelica y genica   *
    ! para los marcadores *
      call diversid(marcadores,n_SNP,n_SNP,divgenmar,divalmar,totalmar)
    call diversid(unobserved,n_neut,n_neut,divgenneut,divalneut,totalneut)

    !   print*,'diversidades marcador'
    !	print*,divalmar,divgenmar

    ! Acumula resultados *
                                  
    ! Calcula parentesco molecular genomico *
      call parmol(fmol,pop,n_pos)

      call promedia_par(fmol,fem,conm)

    !print*,'parentesco genomico'
    !do i=1,npob+1
    !print*,fem(i,:)
    !enddo
    !print*,conm

      call promedia_par(fped,fpm,conped)


    !print*,'parentesco genealogico'
    !do i=1,npob+1
    !print*,fpm(i,:)
    !enddo
    !print*,conped


    ! Calcula parentesco molecular para marcadores *
      call parmol(fmol,marcadores,n_SNP)

      call promedia_par(fmol,fmar,conmar)
      call parmol(fmoln,unobserved,n_neut)

      call promedia_par(fmoln,fmarn,conmarn)


    !print*,'parentesco marcadores'
    !do i=1,npob+1
    !print*,fmar(i,:)
    !enddo
    !print*,conmar

    u=0

    do i=1,sum(ntotal(1,:))
      u=u+fped(i,i)
    enddo
    u=u/sum(ntotal(1,:))
    u=2*u-1
    
    




  results(0,1)=results(0,1)+conped(npob+1)
  results(0,8)=results(0,8)+(conped(npob+1)**2)
  results(0,2)=results(0,2)+fem(npob+1,npob+1)
  results(0,9)=results(0,9)+(fem(npob+1,npob+1)**2)
  results(0,3)=results(0,3)+conm(npob+1)
  results(0,10)=results(0,10)+(conm(npob+1)**2)
  results(0,4)=results(0,4)+homgenom
  results(0,11)=results(0,11)+(homgenom**2)
  results(0,5)=results(0,5)+ibdmar
  results(0,12)=results(0,12)+(ibdmar**2)
  results(0,6)=results(0,6)+dival
  results(0,13)=results(0,13)+(dival**2)
  results(0,7)=results(0,7)+divgen
  results(0,14)=results(0,14)+(divgen**2)
  results(0,15)=results(0,15)+fpm(npob+1,npob+1)
  results(0,16)=results(0,16)+(fpm(npob+1,npob+1)**2)
  results(0,17)=results(0,17)+divalobs
  results(0,18)=results(0,18)+(divalobs**2)
  results(0,19)=results(0,19)+divgenobs
  results(0,20)=results(0,20)+(divgenobs**2)
  results(0,25)=results(0,25)+homar
  results(0,26)=results(0,26)+(homar**2)
  results(0,27)=results(0,27)+fmar(npob+1,npob+1)
  results(0,29)=results(0,29)+(fmar(npob+1,npob+1)**2)
  results(0,28)=results(0,28)+conmar(npob+1)
  results(0,30)=results(0,30)+(conmar(npob+1)**2)
  results(0,34)=results(0,34)+fem(npob+1,npob)
  results(0,37)=results(0,37)+(fem(npob+1,npob)**2)
  results(0,35)=results(0,35)+fpm(npob+1,npob)
  results(0,38)=results(0,38)+(fpm(npob+1,npob)**2)
  results(0,36)=results(0,36)+fmar(npob+1,npob)
  results(0,39)=results(0,39)+(fmar(npob+1,npob)**2)
  results(0,40)=results(0,40)+u
  results(0,41)=results(0,41)+(u**2)

  results2(0,1)=results2(0,1)+conped(npob+1)
  results2(0,8)=results2(0,8)+(conped(npob+1)**2)
  results2(0,2)=results2(0,2)+fem(npob+1,npob+1)
  results2(0,9)=results2(0,9)+(fem(npob+1,npob+1)**2)
  results2(0,3)=results2(0,3)+conm(npob+1)
  results2(0,10)=results2(0,10)+(conm(npob+1)**2)
  results2(0,4)=results2(0,4)+homgenom
  results2(0,11)=results2(0,11)+(homgenom**2)
  results2(0,5)=results2(0,5)+ibdmarn
  results2(0,12)=results2(0,12)+(ibdmarn**2)
  results2(0,17)=results2(0,17)+divalobsn
  results2(0,18)=results2(0,18)+(divalobsn**2)
  results2(0,19)=results2(0,19)+divgenobsn
  results2(0,20)=results2(0,20)+(divgenobsn**2)
  results2(0,25)=results2(0,25)+homarn
  results2(0,26)=results2(0,26)+(homarn**2)
  results2(0,27)=results2(0,27)+fmarn(npob+1,npob+1)
  results2(0,29)=results2(0,29)+(fmarn(npob+1,npob+1)**2)
  results2(0,28)=results2(0,28)+conmarn(npob+1)
  results2(0,30)=results2(0,30)+(conmarn(npob+1)**2)
  results2(0,36)=results2(0,36)+fmarn(npob+1,npob)
  results2(0,39)=results2(0,39)+(fmarn(npob+1,npob)**2)
  results2(0,40)=results2(0,40)+u
  results2(0,41)=results2(0,41)+(u**2)
  
do j=1,npob
  results3(j,0,1)=results3(j,0,1)+homgenoms(j)
  results3(j,0,2)=results3(j,0,2)+ibdmars(j)
  results3(j,0,3)=results3(j,0,3)+homars(j)
  results3(j,0,4)=results3(j,0,4)+homgenomns(j)
  results3(j,0,5)=results3(j,0,5)+ibdmarns(j)
  results3(j,0,6)=results3(j,0,6)+homarns(j)
enddo

!open(unit=432,file='SNP_freq_gen_real.txt')
!    write(432,'(a50,i5)') 'Real population replica number',irep
!open(unit=435,file='Neut_freq_gen_real.txt')
!    write(435,'(a50,i5)') 'Real population replica number',irep
!    
    
if(tfreq.eq.0)      then
    !Calcula frecuencias alelicas  
    allocate(contB(maxval(n_SNP)),pOBSS(n_SNP_total))

    l=0
    contB=0    
    pOBSS=0
    do j=1,ncrom
        do k=1,n_SNP(j)
            l=l+1
            contB(k)=count((marcadores(:,k,j,:))==1)
            pOBSS(l)=dble(contB(k))/dble((nind)*2)
            !write(432,'(i10,f15.10)') l, pOBSS(l)
        end do   
    end do
    
    
    deallocate(contB)
endif

if(tfreq.eq.2)    then  
    !Calcula frecuencias alelicas  
    allocate(contB(maxval(n_SNP)),pOBSSpop(n_SNP_total,npob),pOBSSpop2(n_SNP_total,npob,npob))
    contB=0
    pOBSSpop=0
    kk=1
    o=0    
    
    do i=1,npob
        !write(432,'(a50,i10)') 'subpop', i
        l=0
        o=o+(sum(ntotal(i,:)))
        contB=0    
        do j=1,ncrom
            do k=1,n_SNP(j)
                l=l+1
                contB(k)=count((marcadores(:,k,j,kk:o))==1)
                pOBSSpop(l,i)=dble(contB(k))/dble((sum(ntotal(i,:)))*2)
                !write(432,'(i10,f15.10)') l, pOBSSpop(l,i)
            end do   
        end do
        kk=kk+sum(ntotal(i,:))
    end do
    
    
    
    do i=1,npob
        do j=i,npob
            if(i.ne.j) then
                pOBSSpop2(:,i,j)=(pOBSSpop(:,i)+pOBSSpop(:,j))/2
            else    
                pOBSSpop2(:,i,j)=pOBSSpop(:,i)
            endif
        enddo
    enddo
    
    
  
    
    deallocate(contB,pOBSSpop)
endif  
  

print*,'write freq start'
allocate(contB(maxval(n_SNP)),pOBSSpop(n_SNP_total,npob))
allocate(contN(maxval(n_neut)),pOBSNpop(n_neut_total,npob))
    contB=0
    contN=0
    pOBSNpop=0
    pOBSSpop=0
    kk=1
    o=0 
    do i=1,npob
        !write(432,'(a50,i10)') 'subpop', i
        l=0
        o=o+(sum(ntotal(i,:)))
        contB=0    
        do j=1,ncrom
            !print*, 'ncromsnp',j
            do k=1,n_SNP(j)
                l=l+1
                contB(k)=count((pop(:,SNP_pos(j,k),j,kk:o))==1)
                pOBSSpop(l,i)=dble(contB(k))/dble((sum(ntotal(i,:)))*2)
               ! write(432,'(i10,f15.10)') l, pOBSSpop(l,i)
            end do   
        end do
        kk=kk+sum(ntotal(i,:))
    end do
    
    kk=1
    o=0     
    do i=1,npob
       ! write(435,'(a50,i10)') 'subpop', i
        l=0
        o=o+(sum(ntotal(i,:)))
        contB=0    
        do j=1,ncrom
            !print*, 'ncromneut',j
            do k=1,n_neut(j)
                l=l+1
                contN(k)=count((pop(:,neut_pos(j,k),j,kk:o))==1)
                pOBSNpop(l,i)=dble(contN(k))/dble((sum(ntotal(i,:)))*2)
               ! write(435,'(i10,f15.10)') l, pOBSNpop(l,i)
            end do   
        end do
        kk=kk+sum(ntotal(i,:))
    end do

        
deallocate(contB,pOBSSpop)
        
deallocate(contN,pOBSNpop) 
print*,'write freq end' 
    
! Bucle para generaciones reales*
  do igen=1,ngen
    print*,'Generacion',igen
    
    fgeno=0

    
    
    
 if(tfreq.eq.1)      then
    !Calcula frecuencias alelicas  
    allocate(contB(maxval(n_SNP)),pOBSSpop(n_SNP_total,npob),pOBSSpop2(n_SNP_total,npob,npob))
    pOBSSpop=0
    kk=1
    o=0    
    
    do i=1,npob
        !write(432,'(a50,i10)') 'subpop', i
        l=0
        o=o+(sum(ntotal(i,:)))
        contB=0    
        do j=1,ncrom
            do k=1,n_SNP(j)
                l=l+1
                contB(k)=count((marcadores(:,k,j,kk:o))==1)
                pOBSSpop(l,i)=dble(contB(k))/dble((sum(ntotal(i,:)))*2)
                !write(432,'(i10,f15.10)') l, pOBSSpop(l,i)
            end do   
        end do
        kk=kk+sum(ntotal(i,:))
    end do
    
    
    
    do i=1,npob
        do j=i,npob
            if(i.ne.j) then
                pOBSSpop2(:,i,j)=(pOBSSpop(:,i)+pOBSSpop(:,j))/2
            else    
                pOBSSpop2(:,i,j)=pOBSSpop(:,i)
            endif
        enddo
    enddo
    
    
    deallocate(contB,pOBSSpop)
 endif
  
 print*, 'calculando matriz'
 !!! Calcula matriz de parentesco 
    if(matriz==0) then
                                                    !!!!!!matriz pedigri
        fgeno=fped

            do i=1,nind
                do j=1,nind
                    write(84,*) fped(i,j)
                enddo
            enddo
        
    elseif(matriz==1) then
                                                    !!!!!!matriz similitud
        fgeno=fmol        

    
    elseif(matriz==2) then
                                                    !!!!!!matriz L&H
        if(tfreq.eq.3.or.tfreq.eq.2) then
            call F_and_f_molecularAllo2(nind,npob,ntotal,n_SNP,ncrom,marcadores,pOBSSpop2,fgeno)
        else
            call F_and_f_molecularAllo(nind,n_SNP,ncrom,marcadores,pOBSS,fgeno)         
        endif
        
       
            do i=1,nind
                do j=1,nind
                    write(84,*) fgeno(i,j)
                enddo
            enddo
        
    elseif(matriz==3) then                          !!!!!!matriz VR1
        kkkk=1
        if(tfreq.eq.3.or.tfreq.eq.2) then
            call Gmatrix3(kkkk,npob,ntotal,n_snp,ncrom,nind,marcadores,pOBSSpop2,fgeno)
        else
            call Gmatrix(kkkk,n_SNP,ncrom,nind,marcadores,pOBSS,fgeno)    
        endif
        
               
            do i=1,nind
                do j=1,nind
                    write(84,*) fgeno(i,j)
                enddo
            enddo
            
    elseif(matriz==4) then                          !!!!!!matriz VR2
        kkkk=2
        if(tfreq.eq.3.or.tfreq.eq.2) then
            call Gmatrix3(kkkk,npob,ntotal,n_snp,ncrom,nind,marcadores,pOBSSpop2,fgeno)
        else
            call Gmatrix(kkkk,n_SNP,ncrom,nind,marcadores,pOBSS,fgeno)    
        endif
        
               
            do i=1,nind
                do j=1,nind
                    write(84,*) fgeno(i,j)
                enddo
            enddo
            
            
    elseif(matriz==5) then                          !!!!!!matriz yang
        kkkk=3
        if(tfreq.eq.3.or.tfreq.eq.2) then
            call Gmatrix3(kkkk,npob,ntotal,n_snp,ncrom,nind,marcadores,pOBSSpop2,fgeno)
        else
            call Gmatrix(kkkk,n_SNP,ncrom,nind,marcadores,pOBSS,fgeno)    
        endif
               
            do i=1,nind
                do j=1,nind
                    write(84,*) fgeno(i,j)
                enddo
            enddo
            
    endif
           
if(tfreq.eq.1)      then
    deallocate(pOBSS)
elseif(tfreq.eq.3)    then 
    deallocate(pOBSSpop2)
endif
    
! Calcula los seleccionados con menor probabilidad *
! de perder alelos en la descendencia              *
print*,'annealing'
    call annealing()

! Promedia las contribuciones de cada subpoblacion a las demas *
!write(35,*) 'Generacion',igen
    do i=1,npob
	  pos1=estruc(i,1)
	  ran1=sum(ntotal(i,:))-1
	  do j=1,npob
	    movil(i,j,igen)=movil(i,j,igen)+sum(solopt(pos1:pos1+ran1,j,:))
	  enddo
!write(35,*) movil(i,:,igen),fpm(i,1:npob),conped(i)
	enddo

! Genera genealogia de padres *
    pos=0
    do m=1,npob
      do j=1,2
  	    do i=1,npob
	      pos1=estruc(i,1)
	      ran1=ntotal(i,1)-1
          do k=pos1,pos1+ran1
            if(solopt(k,m,j).gt.0) then
              do l=1,solopt(k,m,j)
                pos=pos+1
!print*,pos,j,k
                geneal(pos,1)=k
              enddo
            endif
          enddo
        enddo
	  enddo
	enddo
          
!print*

! Genera vector de madres *
    pos=0
    do m=1,npob
      do j=1,2
	    do i=1,npob
	      pos1=estruc(i,2)
	      ran1=ntotal(i,2)-1
          do k=pos1,pos1+ran1
            if(solopt(k,m,j).gt.0) then
              do l=1,solopt(k,m,j)
                pos=pos+1
!print*,pos,j,k
                aux(pos)=k
              enddo
            endif
          enddo
        enddo
	  enddo
	enddo

	do i=1,npob
	  do j=1,2
	    pos1=estruc(i,j)
		ran1=ntotal(i,j)-1
		call desordena(aux(pos1:pos1+ran1),ran1+1)
	  enddo
	enddo

! Construye genealogia via madres *
    do i=1,nind
      geneal(i,2)=aux(i)
    enddo



print*, 'genera hijos generacion', igen
! Genera nueva poblacion *
    do i=1,nind
      do j=1,ncrom
        do k=1,2
! Decide el numero de sobrecruzamientos y su posicion *
! y los ordena                                        *
 ! For small chromosomes (less than 1 Morgan) !
            if(long.lt.1) then
	            call mother(u)
            
	            if(u.lt.long) then
	                long=1
	            else
		            long=0
                endif
            else
                long=long0
            
            endif
            
              nsob=0
              do n=1,anint(long)
    
                call mother(u)
                !print*,u,'Entra en Poisson'
                call poisson(u*.9999,l)
                nsob=nsob+l
              enddo
                
              deallocate (cross)
              allocate (cross(nsob))
              cross=0
              do n=1,nsob
                call mother(u)
                cross(n)=min(int(u*genom_length(j)+1),genom_length(j))
                
		        cross(n)=anint(cross(n)/interval(j))
              enddo

              if(nsob.gt.1) call ordena(cross,nsob)

    ! Copia las partes de los homologos correspondientes *
              call mother(u)
              pad=min(int(u*2+1),2)
              sobant=0
              pos=1
              do n=1,nsob
                popt1(k,sobant+1:cross(n),j,i)=pop(pad,sobant+1:cross(n),j,geneal(i,k))
                sobant=cross(n)
                do l=pos,n_SNP(j)
                  if(SNP_pos(j,l).ge.cross(n)) exit
	              marct1(k,l,j,i)=marcadores(pad,l,j,geneal(i,k))
	              pos=pos+1
	            enddo
                select case(pad)
                  case(1)
                    pad=2
                  case(2)
                    pad=1
                end select
              enddo
              popt1(k,sobant+1:n_pos(j),j,i)=pop(pad,sobant+1:n_pos(j),j,geneal(i,k))
	            if(pos.le.n_SNP(j)) then
                   marct1(k,pos:n_SNP(j),j,i)=marcadores(pad,pos:n_SNP(j),j,geneal(i,k))
	            endif
            enddo
          enddo
        enddo
    pop=popt1
    popt1=0
	marcadores=marct1
	marct1=0
    unobserved=0
    
    do i=1,nind
          do j=1,ncrom
              do l=1,n_neut(j)
                      unobserved(:,l,j,i)=pop(:,neut_pos(j,l),j,i)  
              enddo
          enddo
      enddo  


! Calcula parentesco genealogico *
    do i=1,nind
      fpedt1(i,i)=.5*(1+fped(geneal(i,1),geneal(i,2)))
      do j=i+1,nind
        fpedt1(i,j)=.25*(fped(geneal(i,1),geneal(j,1))+&
                         fped(geneal(i,1),geneal(j,2))+&
                         fped(geneal(i,2),geneal(j,1))+&
                         fped(geneal(i,2),geneal(j,2)))
        fpedt1(j,i)=fpedt1(i,j)
      enddo
    enddo
    fped=fpedt1

    call promedia_par(fped,fpm,conped)



!print*,'parentesco genealogico'
!do i=1,npob+1
!print*,fpm(i,:)
!enddo
!print*,conped


! Calcula heterocigosidad real *
          call hetero(pop,n_pos,n_SNP,SNP_pos,marcadores,homgenom,ibdmar,homar)
    call hetero(pop,n_pos,n_neut,neut_pos,unobserved,homgenomn,ibdmarn,homarn)
          call hetero2(pop,n_pos,n_SNP,SNP_pos,marcadores,homgenoms,ibdmars,homars)
    call hetero2(pop,n_pos,n_neut,neut_pos,unobserved,homgenomns,ibdmarns,homarns)

!     print*,'homocigosidad genomica,pos. marcadores, marcadores'
!	print*,homgenom,ibdmar,homar

! Calcula diversidad alelica y genica *
          call diversid(pop,n_pos,n_pos,divgen,dival,total)

!      print*,'diversidades genomicas'
!	print*,dival,divgen

! Calcula diversidad alelica y genica   *
! para las posiciones de los marcadores *
          call diversid(pop,n_SNP,n_pos,divgenobs,divalobs,SNP_pos)
    call diversid(pop,n_neut,n_pos,divgenobsn,divalobsn,neut_pos)
                 
!      print*,'diversidades pos. marcador'
!	print*,divalobs,divgenobs

! Calcula diversidad alelica y genica   *
! para los marcadores *
          call diversid(marcadores,n_SNP,n_SNP,divgenmar,divalmar,totalmar)
    call diversid(unobserved,n_neut,n_neut,divgenneut,divalneut,totalneut)
                 
!      print*,'diversidades marcador'
!	print*,divalmar,divgenmar

! Acumula resultados *

! Calcula parentesco molecular genomico *
    call parmol(fmol,pop,n_pos)

    call promedia_par(fmol,fem,conm)

!print*,'parentesco genomico'
!do i=1,npob+1
!print*,fem(i,:)
!enddo
!print*,conm



! Calcula parentesco molecular para marcadores *
    call parmol(fmol,marcadores,n_SNP)

    call promedia_par(fmol,fmar,conmar)
      call parmol(fmoln,unobserved,n_neut)

      call promedia_par(fmoln,fmarn,conmarn)


!print*,'parentesco marcadores'
!do i=1,npob+1
!print*,fmar(i,:)
!enddo
!print*,conmar

u=0

do i=1,sum(ntotal(1,:))
  u=u+fped(i,i)
enddo
u=u/sum(ntotal(1,:))
u=2*u-1
print*, 'results'
  results(igen,1)=results(igen,1)+conped(npob+1)
  results(igen,8)=results(igen,8)+(conped(npob+1)**2)
  results(igen,2)=results(igen,2)+fem(npob+1,npob+1)
  results(igen,9)=results(igen,9)+(fem(npob+1,npob+1)**2)
  results(igen,3)=results(igen,3)+conm(npob+1)
  results(igen,10)=results(igen,10)+(conm(npob+1)**2)
  results(igen,4)=results(igen,4)+homgenom
  results(igen,11)=results(igen,11)+(homgenom**2)
  results(igen,5)=results(igen,5)+ibdmar
  results(igen,12)=results(igen,12)+(ibdmar**2)
  results(igen,6)=results(igen,6)+dival
  results(igen,13)=results(igen,13)+(dival**2)
  results(igen,7)=results(igen,7)+divgen
  results(igen,14)=results(igen,14)+(divgen**2)
  results(igen,15)=results(igen,15)+fpm(npob+1,npob+1)
  results(igen,16)=results(igen,16)+(fpm(npob+1,npob+1)**2)
  results(igen,17)=results(igen,17)+divalobs
  results(igen,18)=results(igen,18)+(divalobs**2)
  results(igen,19)=results(igen,19)+divgenobs
  results(igen,20)=results(igen,20)+(divgenobs**2)
  results(igen,25)=results(igen,25)+homar
  results(igen,26)=results(igen,26)+(homar**2)
  results(igen,27)=results(igen,27)+fmar(npob+1,npob+1)
  results(igen,29)=results(igen,29)+(fmar(npob+1,npob+1)**2)
  results(igen,28)=results(igen,28)+conmar(npob+1)
  results(igen,30)=results(igen,30)+(conmar(npob+1)**2)
  results(igen,31)=results(igen,31)+particion(1,2)
  results(igen,32)=results(igen,32)+particion(2,2)
  results(igen,33)=results(igen,33)+particion(3,2)
  results(igen,34)=results(igen,34)+fem(npob+1,npob)
  results(igen,37)=results(igen,37)+(fem(npob+1,npob)**2)
  results(igen,35)=results(igen,35)+fpm(npob+1,npob)
  results(igen,38)=results(igen,38)+(fpm(npob+1,npob)**2)
  results(igen,36)=results(igen,36)+fmar(npob+1,npob)
  results(igen,39)=results(igen,39)+(fmar(npob+1,npob)**2)
  results(igen,40)=results(igen,40)+u
  results(igen,41)=results(igen,41)+(u**2)
 

  results2(igen,1)=results2(igen,1)+conped(npob+1)
  results2(igen,8)=results2(igen,8)+(conped(npob+1)**2)
  results2(igen,2)=results2(igen,2)+fem(npob+1,npob+1)
  results2(igen,9)=results2(igen,9)+(fem(npob+1,npob+1)**2)
  results2(igen,3)=results2(igen,3)+conm(npob+1)
  results2(igen,10)=results2(igen,10)+(conm(npob+1)**2)
  results2(igen,4)=results2(igen,4)+homgenom
  results2(igen,11)=results2(igen,11)+(homgenom**2)
  results2(igen,5)=results2(igen,5)+ibdmarn
  results2(igen,12)=results2(igen,12)+(ibdmarn**2)
  results2(igen,17)=results2(igen,17)+divalobsn
  results2(igen,18)=results2(igen,18)+(divalobsn**2)
  results2(igen,19)=results2(igen,19)+divgenobsn
  results2(igen,20)=results2(igen,20)+(divgenobsn**2)
  results2(igen,25)=results2(igen,25)+homarn
  results2(igen,26)=results2(igen,26)+(homarn**2)
  results2(igen,27)=results2(igen,27)+fmarn(npob+1,npob+1)
  results2(igen,29)=results2(igen,29)+(fmarn(npob+1,npob+1)**2)
  results2(igen,28)=results2(igen,28)+conmarn(npob+1)
  results2(igen,30)=results2(igen,30)+(conmarn(npob+1)**2)
  results2(igen,36)=results2(igen,36)+fmarn(npob+1,npob)
  results2(igen,39)=results2(igen,39)+(fmarn(npob+1,npob)**2)
  results2(igen,40)=results2(igen,40)+u
  results2(igen,41)=results2(igen,41)+(u**2)
  
do j=1,npob
  results3(j,igen,1)=results3(j,igen,1)+homgenoms(j)
  results3(j,igen,2)=results3(j,igen,2)+ibdmars(j)
  results3(j,igen,3)=results3(j,igen,3)+homars(j)
  results3(j,igen,4)=results3(j,igen,4)+homgenomns(j)
  results3(j,igen,5)=results3(j,igen,5)+ibdmarns(j)
  results3(j,igen,6)=results3(j,igen,6)+homarns(j)
enddo
  
  
  print*,'write freq start'
allocate(contB1(maxval(n_SNP)),pOBSSpop1(n_SNP_total,npob))
allocate(contN1(maxval(n_neut)),pOBSNpop1(n_neut_total,npob))
    pOBSNpop1=0
    pOBSSpop1=0
    contB1=0
    contN1=0
    kk=1
    o=0 
  !  write(432,'(a50,i10)') 'Generacion', igen
    do i=1,npob
        !write(432,'(a50,i10)') 'subpop', i
        l=0
        o=o+(sum(ntotal(i,:)))
        contB1=0    
        do j=1,ncrom
            do k=1,n_SNP(j)
                l=l+1
                contB1(k)=count((pop(:,SNP_pos(j,k),j,kk:o))==1)
                pOBSSpop1(l,i)=dble(contB1(k))/dble((sum(ntotal(i,:)))*2)
              !  write(432,'(i10,f15.10)') l, pOBSSpop1(l,i)
            end do   
        end do
        kk=kk+sum(ntotal(i,:))
    end do
    
    kk=1
    o=0 
   ! write(435,'(a50,i10)') 'Generacion', igen
    do i=1,npob
       ! write(435,'(a50,i10)') 'subpop', i
        l=0
        o=o+(sum(ntotal(i,:)))
        contB1=0    
        do j=1,ncrom
            do k=1,n_neut(j)
                l=l+1
                contN1(k)=count((pop(:,neut_pos(j,k),j,kk:o))==1)
                pOBSNpop1(l,i)=dble(contN1(k))/dble((sum(ntotal(i,:)))*2)
                !write(435,'(i10,f15.10)') l, pOBSNpop1(l,i)
            end do   
        end do
        kk=kk+sum(ntotal(i,:))
    end do

        
deallocate(contB1,pOBSSpop1)
        
deallocate(contN1,pOBSNpop1) 
print*,'write freq end'
  
! Cierra bucle de generaciones reales*
  enddo
print*, 'termina bucle generaciones de replica',irep
print*, ' Genera salida de resultados '

  call cpu_time(time3)


deallocate(SNP_pos,genom_map,neut_pos,genom_map_neut)
deallocate (totalmar,total,totalneut)
deallocate (pop,popt1,marcadores,marct1,unobserved)

 open(5,FILE="meta_gen.sal")
  write(5,*) 'Metapoblaciones con parentesco genealogico'
  write(5,*) 'N. individuos',nind
  write(5,*) 'N. poblaciones',npob
  write(5,*) 'Tamaño poblaciones',ntotal(:,1)
  write(5,*) '                  ',ntotal(:,2)
  write(5,*) 'Peso diversidad dentro',lambda,'N. maximo de migraciones',nmigra/2
  write(5,*) 'N. cromosomas', ncrom,' Longitud',long
  write(5,*) 'N. loci por cromosoma',n_pos,' con alelos'
  write(5,*) 'N. de marcadores genotipados',n_SNP,' con alelos'
  write(5,*) 'N. generaciones sin control',ngs
  write(5,*) 'N. generaciones',ngen,' N. repeticiones',nrep
  write(5,*) 'Matriz usada',matriz
  write(5,*) '*************************************'
  write(5,*) 'Repeticion',irep,'Tiempo',time3-time1,time3-time2
  time2=time3
  write(5,*) 
        write(5,'(a2,7a10)') '  ','conped','fbgenom','conm','fbped',&
     'entre','dentro','migra'
        do i=0,ngen
          write(5,'(i2,7f10.6)') i,(results(i,j)/irep,j=1,3),&
                           results(i,15)/irep,(results(i,j)/irep,j=31,33)
          write(5,'(a2,4f10.6)') '  ',((results(i,j+7)/irep)-&
                           (results(i,j)/irep)**2,j=1,3),&
     (results(i,16)/irep)-(results(i,15)/irep)**2
        enddo
        write(5,*)
        write(5,'(a2,9a10)') '  ','homgeno',' ibdmar','dival',&
     'divgen','divalobs','divgenobs','homar','fbmar','conmar'
        do i=0,ngen
          write(5,'(i2,9f10.6)') i,(results(i,j)/irep,j=4,7),&
                              (results(i,j)/irep,j=17,19,2),&
                              results(i,25)/irep,&
                              (results(i,j)/irep,j=27,28)
          write(5,'(a2,9f10.6)') '  ',((results(i,j+7)/irep)-&
                           (results(i,j)/irep)**2,j=4,7),&
                           ((results(i,j+1)/irep)-&
                           (results(i,j)/irep)**2,j=17,19,2),&
                           ((results(i,26)/irep)-&
                           (results(i,25)/irep)**2),&
                           ((results(i,j+2)/irep)-&
                           (results(i,j)/irep)**2,j=27,28)
        enddo
        write(5,*) 
        write(5,'(a2,3a10)') '  ','fwgenom','fwped','fwmar'
        do i=0,ngen
          write(5,'(i2,7f10.6)') i,(results(i,j)/irep,j=34,36)
          write(5,'(a2,4f10.6)') '  ',((results(i,j+3)/irep)-&
                           (results(i,j)/irep)**2,j=34,36)
        enddo
        write(5,*) 
        write(5,'(a2,1a10)') '  ','fped_1'
        do i=0,ngen
          write(5,'(i2,f10.6)') i,results(i,40)/irep
          write(5,'(a2,f10.6)') '  ',(results(i,41)/irep)-&
                           (results(i,40)/irep)**2
        enddo
        
        
         write(5,*) 
        write(5,'(a2,7a10)') '  ','conpedn','fbgenomn','conmn','fbped',&
     'entre','dentro','migra'
        do i=0,ngen
          write(5,'(i2,7f10.6)') i,(results2(i,j)/irep,j=1,3),&
                           results2(i,15)/irep,(results2(i,j)/irep,j=31,33)
          write(5,'(a2,4f10.6)') '  ',((results2(i,j+7)/irep)-&
                           (results2(i,j)/irep)**2,j=1,3),&
     (results2(i,16)/irep)-(results2(i,15)/irep)**2
        enddo
        write(5,*)
        write(5,'(a2,9a12)') '  ','homgenon',' ibdmarn','divalkk',&
     'divgenkk','divalobs','divgenobsn','homarn','fbmarn','conmarn'
        do i=0,ngen
          write(5,'(i2,9f10.6)') i,(results2(i,j)/irep,j=4,7),&
                              (results2(i,j)/irep,j=17,19,2),&
                              results2(i,25)/irep,&
                              (results2(i,j)/irep,j=27,28)
          write(5,'(a2,9f10.6)') '  ',((results2(i,j+7)/irep)-&
                           (results2(i,j)/irep)**2,j=4,7),&
                           ((results2(i,j+1)/irep)-&
                           (results2(i,j)/irep)**2,j=17,19,2),&
                           ((results2(i,26)/irep)-&
                           (results2(i,25)/irep)**2),&
                           ((results2(i,j+2)/irep)-&
                           (results2(i,j)/irep)**2,j=27,28)
        enddo
        write(5,*) 
        write(5,'(a2,3a10)') '  ','fwgenom','fwped','fwmar'
        do i=0,ngen
          write(5,'(i2,7f10.6)') i,(results2(i,j)/irep,j=34,36)
          write(5,'(a2,4f10.6)') '  ',((results2(i,j+3)/irep)-&
                           (results2(i,j)/irep)**2,j=34,36)
        enddo
        write(5,*) 
        write(5,'(a2,1a10)') '  ','fped_1'
        do i=0,ngen
          write(5,'(i2,f10.6)') i,results2(i,40)/irep
          write(5,'(a2,f10.6)') '  ',(results2(i,41)/irep)-&
                           (results2(i,40)/irep)**2
        enddo
        
        
        do i=0,ngen
            do k=1,npob
                write(5,'(a2,4a10)') '  ','gen','subpop','homgenos',' ibdmars',&
                    'homars','homgenons',' ibdmarns','homarns'
                write(5,'(2i2,6f10.6)') i,k,(results3(k,i,j)/irep,j=1,6)
            enddo
        enddo      
       
        close(5)
        
        
    open(35,FILE="meta_contribuciones.txt")
    do i=1,ngen
      write(35,*) 'Generación',i
      write(35,*) (j,j=1,npob)
      do j=1,npob
        write(35,'(i2,10f6.2)') j,dble(movil(j,:,i))/(irep*2)
      enddo
    enddo
    
    
    close(35)
    
if(tfreq.eq.0)      then
    deallocate(pOBSS)
elseif(tfreq.eq.2)    then 
    deallocate(pOBSSpop2)
endif
! Cierra bucle de repeticiones *
enddo



print*,'PROGRAMA FINALIZADO'


      end program

subroutine annealing()
! Minimiza la funcion usando simulated annealing *
use comun
implicit none
integer sol,sola,ch,pos,anc,i,j,l,n,m,r,s1,s2,p1,p2,ivig,niv,rep,camb,pos1,ran1
double precision kt,k,t,delta,omega,u,eneac,eneal,eneopt
parameter(t=.01,k=.9)
dimension sol(nind,npob,2),sola(nind,npob,2)

! Solucion inicial al azar *      
sol=0
do i=1,npob
  do j=1,2
    do l=1,ntotal(i,j)
	  call mother(u)
	  n=min(int(u*nind)+1,nind)
	  sol(n,i,j)=sol(n,i,j)+1
	  do m=2,npob
	    p1=m-1
		if(n.lt.estruc(m,1)) exit
		p1=npob
	  enddo
	  s1=2
	  if(n.ge.estruc(p1,2)) s1=1
	  pos1=estruc(p1,s1)
	  ran1=ntotal(p1,s1)
	  call mother(u)
	  m=min(int(u*ran1)+pos1,pos1+ran1-1)
	  sol(m,i,j)=sol(m,i,j)+1
    enddo
  enddo
enddo
kt=t/k      

! Calcula valor de estado inicial *      
call energia(sol,eneac)
!print*,eneac
solopt=sol
eneopt=eneac

ivig=5000

! Bucle de niveles *     
do niv=1,50
  camb=1+(10*ivig/5000)
  kt=kt*k
  ivig=0

! Bucle de repeticiones *      
  do rep=1,5000
! Genera solucion alternativa *
! Copia la solucion actual en la alternativa *
    sola=sol
	do r=1,camb
! Decide el individuo (y su pareja) que pierden un hijo *
! Primero el individuo !
      individuo:do
	    call mother(u)
		n=min(int(u*nind)+1,nind)
		if(sum(sola(n,:,:)).eq.0) cycle individuo
		poblacion:do
		  call mother(u)
		  p1=min(int(u*npob)+1,npob)
		  if(sum(sola(n,p1,:)).eq.0) cycle poblacion
          call mother(u)
          s1=min(int(u*2)+1,2)
		  if(sola(n,p1,s1).eq.0) then
		    select case(s1)
		      case(1)
		        s1=2
			  case(2)
			    s1=1
			endselect
		  endif
		  sola(n,p1,s1)=sola(n,p1,s1)-1
		  exit individuo
	    enddo poblacion
	  enddo individuo

! y mira de que poblacion y sexo es !

	  do m=2,npob
	    p2=m-1
		if(n.lt.estruc(m,1)) exit
		p2=npob
	  enddo
	  s2=2
	  if(n.ge.estruc(p2,2)) s2=1

! Ahora elige una pareja de la misma poblacion !
      pos1=estruc(p2,s2)
	  ran1=ntotal(p2,s2)
	  do
	    call mother(u)
		m=min(int(u*ran1)+pos1,pos1+ran1-1)
        if(sola(m,p1,s1).eq.0) cycle
  	    sola(m,p1,s1)=sola(m,p1,s1)-1
		exit
	  enddo

! Decide el individuo (y su pareja) que ganan un hijo *
! Primero el individuo !
	  call mother(u)
	  p2=min(int(u*npob)+1,npob)
      pos1=estruc(p2,s2)
	  ran1=ntotal(p2,s2)
      do
	    call mother(u)
		l=min(int(u*ran1)+pos1,pos1+ran1-1)
        if(m.eq.l) cycle
  	    sola(l,p1,s1)=sola(l,p1,s1)+1
		exit
	  enddo
	    
! Ahora elige una pareja de la misma poblacion !
      select case(s2)
	    case(1)
		  s2=2
		case(2)
		  s2=1
	  endselect
	  pos1=estruc(p2,s2)
	  ran1=ntotal(p2,s2)
      call mother(u)
      n=min(int(u*ran1)+pos1,pos1+ran1-1)
  	  sola(n,p1,s1)=sola(n,p1,s1)+1
	enddo

! Calcula valor de estado alternativo
    call energia(sola,eneal)
!print*,eneal,eneopt

! Anneling *
    ch=0
    delta=amax1(eneal-eneac,0.d0)
    omega=exp(-delta/kt)
    if(omega.ge.1) then
      ch=1
      if(eneal.lt.eneopt) then
        eneopt=eneal
        solopt=sola
		particion(:,2)=particion(:,1)
      endif
    else
      call mother(u)
      if(u.lt.omega) ch=1
    endif
    if(ch.eq.1) then
      sol=sola
	  eneac=eneal
      ivig=ivig+1
    endif
  enddo
  !print*,niv,ivig,1+(10*ivig/5000),eneac,eneopt
  if(ivig.eq.0) then
    return
  endif
enddo
!write(48,*) 'Bucles agotados'
end subroutine annealing


subroutine energia(sol,ene)
! Calcula valor de la funcion objetivo *
use comun
use decision
implicit none
integer sol,i,j,k,l,s1,s2,pos,ran,migra,penaliza
double precision ene,dentro,entre
dimension sol(nind,npob,2)

ene=0.
penaliza=0
dentro=0
entre=0
migra=sum(sol)
do k=1,npob
  do s1=1,2
    pos=estruc(k,s1)
    ran=ntotal(k,s1)-1
    migra=migra-sum(sol(pos:pos+ran,k,:))
  enddo
  do s1=1,2
    do s2=1,2
	  do i=1,nind
	    dentro=dentro+(sol(i,k,s1)*sol(i,k,s2))*fgeno(i,i)/(4*ntotal(k,s1)*ntotal(k,s2))
	    do j=i+1,nind
	      dentro=dentro+2*sol(i,k,s1)*sol(j,k,s2)*fgeno(i,j)/(4*ntotal(k,s1)*ntotal(k,s2))
	    enddo
	  enddo
	enddo
  enddo
  do l=1,npob
    if(l.eq.k) cycle
    do s1=1,2
      do s2=1,2
	    do i=1,nind
	      entre=entre+(sol(i,k,s1)*sol(i,l,s2))*fgeno(i,i)/(4*ntotal(k,s1)*ntotal(l,s2))
	      do j=i+1,nind
	        entre=entre+2*sol(i,k,s1)*sol(j,l,s2)*fgeno(i,j)/(4*ntotal(k,s1)*ntotal(l,s2))
	      enddo
	    enddo
      enddo
    enddo
  enddo
enddo

particion(1,1)=entre
particion(2,1)=dentro
particion(3,1)=migra

if(migra.gt.nmigra) penaliza=(migra-nmigra)*10000
ene=entre+lambda*dentro+penaliza
end subroutine energia

subroutine parmol(ft,pop,nloci)
! Calcula la matriz de parentescos moleculares (multialelicos) *
use comun
implicit none
integer*1  pop
integer i,j,k1,k2,l,n,iale,icont,nloci
double precision ft
dimension nloci(ncrom)
dimension ft(nind,nind),pop(2,maxval(nloci),ncrom,nind)
dimension iale(4)
ft=0
do i=1,nind
  do j=i,nind
    icont=0
    do l=1,ncrom
      do n=1,nloci(l)
        do k1=1,2
          do k2=1,2
            if(pop(k1,n,l,i).eq.pop(k2,n,l,j)) icont=icont+1
          enddo
        enddo
      enddo
    enddo
    ft(j,i)=dble(icont)/(4*sum(nloci))
    ft(i,j)=ft(j,i)
  enddo
enddo
end subroutine


! Calcula diversidad genica (1 - sum(q^2))    *
! y alelica (numero de alelos supervivientes) *
subroutine diversid(pop,nlusa,nloci,divgen,dival,mpos)
use comun
implicit none
integer*1  pop
integer l,i,j,m,n,nale,iale,k,nalel,mpos,nlusa,nloci,pos,ran
double precision divgen,dival,sumq2
dimension nloci(ncrom)
dimension nlusa(ncrom)
dimension pop(2,maxval(nloci),ncrom,nind),nale(2,2)
dimension mpos(ncrom,maxval(nlusa))

sumq2=0
dival=0
do j=1,ncrom
  do l=1,nlusa(j)
    nale=0
    do i=1,npob
	  do k=1,2
	    pos=estruc(i,k)
		ran=ntotal(i,k)-1
	    do n=pos,pos+ran
          do m=1,2
            iale=pop(m,mpos(j,l),j,n)+1
            nale(iale,k)=nale(iale,k)+1
          enddo
        enddo
      enddo
    enddo
    do k=1,2
      sumq2=sumq2+(dble(nale(k,1))/(4*nm)+dble(nale(k,2))/(4*nh))**2
      if(sum(nale(k,:)).ne.0) dival=dival+1
    enddo
  enddo
enddo
divgen=1-(sumq2/(sum(nlusa)))
dival=dival/(2*sum(nlusa))
    end subroutine

    
! Calcula heterocigosidad real *
subroutine hetero(pop,nloci,nlusa,mpos,marcadores,homgenom,ibdmar,homar)
use comun
implicit none
integer*1  pop,marcadores
integer nloci,mpos,i,j,k,l,nlusa
double precision homgenom,ibdmar,homar
dimension nloci(ncrom)
dimension nlusa(ncrom)
dimension pop(2,maxval(nloci),ncrom,nind),mpos(ncrom,maxval(nlusa)),marcadores(2,maxval(nlusa),ncrom,nind)

ibdmar=0
homgenom=0
homar=0

do i=1,nind
  do j=1,ncrom
    do l=1,nlusa(j)
      if(pop(1,mpos(j,l),j,i).ne.pop(2,mpos(j,l),j,i)) ibdmar=ibdmar+1
      if(marcadores(1,l,j,i).ne.marcadores(2,l,j,i)) homar=homar+1
    enddo
  enddo
enddo

do i=1,nind
  do j=1,ncrom
    do l=1,nloci(j)
      if(pop(1,l,j,i).ne.pop(2,l,j,i)) homgenom=homgenom+1
    enddo
  enddo
enddo

ibdmar=1-(ibdmar/(sum(nlusa)*nind))
homar=1-(homar/(sum(nlusa)*nind))
homgenom=1-(homgenom/(sum(nloci)*nind))
    end subroutine


       
! Calcula heterocigosidad real por subpop*
subroutine hetero2(pop,nloci,nlusa,mpos,marcadores,homgenom,ibdmar,homar)
use comun
implicit none
integer*1  pop,marcadores
integer nloci,mpos,i,j,kk,o,k,l,nlusa
double precision homgenom,ibdmar,homar
dimension nloci(ncrom)
dimension nlusa(ncrom)
dimension pop(2,maxval(nloci),ncrom,nind),mpos(ncrom,maxval(nlusa)),marcadores(2,maxval(nlusa),ncrom,nind)
dimension homgenom(npob),ibdmar(npob),homar(npob)
ibdmar=0
homgenom=0
homar=0

kk=1
o=0 
do k=1,npob
    o=o+(sum(ntotal(k,:)))
do i=1,nind
    if(i.ge.kk.and.i.le.o)then
  do j=1,ncrom
    do l=1,nlusa(j)
      if(pop(1,mpos(j,l),j,i).ne.pop(2,mpos(j,l),j,i)) ibdmar(k)=ibdmar(k)+1
      if(marcadores(1,l,j,i).ne.marcadores(2,l,j,i)) homar(k)=homar(k)+1
    enddo
  enddo
    endif
enddo
ibdmar(k)=1-(ibdmar(k)/(sum(nlusa)*sum(ntotal(k,:))))
homar(k)=1-(homar(k)/(sum(nlusa)*sum(ntotal(k,:))))
    kk=kk+sum(ntotal(k,:))
enddo


kk=1
o=0
do k=1,npob
    o=o+(sum(ntotal(k,:)))
do i=1,nind
    if(i.ge.kk.and.i.le.o)then
  do j=1,ncrom
    do l=1,nloci(j)
      if(pop(1,l,j,i).ne.pop(2,l,j,i)) homgenom(k)=homgenom(k)+1
    enddo
  enddo
    endif
enddo
homgenom(k)=1-(homgenom(k)/(sum(nloci)*sum(ntotal(k,:))))
    kk=kk+sum(ntotal(k,:))
enddo

end subroutine

 
    
    
    ! Aleatoriza un vector *
subroutine desordena(ix,nx)
double precision u
dimension ix(nx)
do i=1,nx
  call mother(u)
  n=u*nx*.999999+1
  ic=ix(i)
  ix(i)=ix(n)
  ix(n)=ic
enddo
end

! Ordena un vector por los valores que contiene *
subroutine ordena(ix,nx)
dimension ix(nx)
do i=1,nx-1
  do j=i+1,nx
    if(ix(j).lt.ix(i)) then
      ic=ix(i)
      ix(i)=ix(j)
      ix(j)=ic
    endif
  enddo
enddo
end

subroutine promedia_par(f,fm,conm)
!Promedia el parentesco entre y dentro de poblaciones !
use comun
implicit none
integer pos1,pos2,ran1,ran2,i,j,k,k1,k2,l
double precision f,fm,conm
dimension f(nind,nind),fm(npob+1,npob+1),conm(npob+1)
fm=0
conm=0

do i=1,npob
  do k=1,2
    pos1=estruc(i,k)
    ran1=ntotal(i,k)-1
    fm(i,i)=fm(i,i)+sum(f(pos1:pos1+ran1,pos1:pos1+ran1))/(ntotal(i,k)**2)
    do l=pos1,pos1+ran1
      conm(i)=conm(i)+(2*f(l,l)-1)/ntotal(i,k)
    enddo
  enddo
  pos1=estruc(i,1)
  ran1=ntotal(i,1)-1
  pos2=estruc(i,2)
  ran2=ntotal(i,2)-1
  fm(i,i)=fm(i,i)+2*sum(f(pos1:pos1+ran1,pos2:pos2+ran2))/(ntotal(i,1)*ntotal(i,2))
  fm(i,i)=fm(i,i)/4
  conm(i)=conm(i)/2
  do j=i+1,npob
    do k1=1,2
      pos1=estruc(i,k1)
      ran1=ntotal(i,k1)-1
      do k2=1,2
        pos2=estruc(j,k2)
        ran2=ntotal(j,k2)-1
        fm(i,j)=fm(i,j)+sum(f(pos1:pos1+ran1,pos2:pos2+ran2))/(ntotal(i,k1)*ntotal(j,k2))
      enddo
    enddo
    fm(i,j)=fm(i,j)/4
    fm(j,i)=fm(i,j)
  enddo
enddo
conm(npob+1)=dot_product(conm(1:npob),sum(ntotal,2))/nind
do i=1,npob
  fm(npob+1,npob)=fm(npob+1,npob)+fm(i,i)*sum(ntotal(i,:))**2
  do j=i+1,npob
    fm(npob+1,npob+1)=fm(npob+1,npob+1)+2*fm(i,j)*sum(ntotal(i,:))*sum(ntotal(j,:))
  enddo
enddo
fm(npob+1,npob)=fm(npob+1,npob)/sum((sum(ntotal,2)**2))
fm(npob+1,npob+1)=fm(npob+1,npob+1)/((nind**2)-sum((sum(ntotal,2)**2)))
end subroutine

! Rutina para numeros aleatorios uniformes (0 - 1) *
      subroutine setup(iseed)
      Common /seeds/xseed(8),yseed(8),i16,i32
      integer iseed,xseed,yseed,carry
      i16=65535
      i32=2147483647
      n=iseed
      do 10 i=1,8
         k=mod(n,1000)
         carry=mod(n/1000,1000)
         n=672*k+carry
         xseed(i)=mod(n,i16)
         n=672*carry+k
         yseed(i)=mod(n,i16)
  10  continue
      end

      subroutine mother(u)
      Common /seeds/x(8),y(8),i16,i32
      integer x,y,carx,cary
      double precision u,w,t
      k=mod(x(8),21475)*100000+y(8)+mod(x(8),18112)
       if(k.lt.0)k=-k
      k=mod(k,i32)
      L=mod(y(8),21475)*100000+x(8)+mod(y(8),18112)
       if(L.lt.0)L=-L
      L=mod(L,i32) 
      t=i32
      carx=mod(L,i16)
      cary=mod(k,i16)
       n=x(1)*12013 + x(2)*1066 + x(3)*1215 + x(4)*1492&
        + x(5)*1776 + x(6)*1812 + x(7)*1860 + x(8)*1941 + carx
       m=y(1)*9272 + y(2)*7777 + y(3)*6666 + y(4)*5555&
        + y(5)*4444 + y(6)*3333 + y(7)*2222 + y(8)*1111 + cary
      do 10 i=1,7
       x(i)=x(i+1)
  10   y(i)=y(i+1)

      x(8)=mod(n,i16)
      y(8)=mod(m,i16)
      k=mod(n,21475)*100000+mod(m,83648)
       if(k.lt.0)k=-k
      k=mod(k,i32)
      w=k
       u=w/t
      end


!  dirichlet.f90  17 June 2000
!
function alngam ( x, ifault )
!
!*******************************************************************************
!
!! ALNGAM computes the logarithm of the Gamma function.
!
!
!  Reference:
!
!    Allan Macleod,
!    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
!    Algorithm AS 245,
!    Applied Statistics,
!    Volume 38, Number 2, pages 397-402, 1989.
!
!  Modified:
!
!    29 November 1999
!
!  Parameters:
!
!    Input, real X, the argument of the Gamma function.
!
!    Output, integer IFAULT, error flag.
!    0, no error occurred.
!    1, X is less than or equal to 0.
!    2, X is too big.
!
!    Output, real ALNGAM, the logarithm of the gamma function of X.
!
  double precision, parameter :: alr2pi = 0.918938533204673
  double precision, parameter :: xlge = 5.10E+6
  double precision, parameter :: xlgst = 1.0E+30
!
  double precision alngam
  integer ifault
  double precision, parameter, dimension ( 9 ) :: r1 = (/ &
     -2.66685511495, &
     -24.4387534237, &
     -21.9698958928, &
      11.1667541262, &
       3.13060547623, &
       0.607771387771, &
      11.9400905721, &
      31.4690115749, &
      15.2346874070 /)
  double precision, parameter, dimension ( 9 ) :: r2 = (/ &
    -78.3359299449, &
    -142.046296688, &
     137.519416416, &
     78.6994924154, &
     4.16438922228, &
     47.0668766060, &
     313.399215894, &
     263.505074721, &
     43.3400022514 /)
  double precision, parameter, dimension ( 9 ) :: r3 = (/ &
    -2.12159572323E5, &
     2.30661510616E5, &
     2.74647644705E4, &
    -4.02621119975E4, &
    -2.29660729780E3, &
    -1.16328495004E5, &
    -1.46025937511E5, &
    -2.42357409629E4, &
    -5.70691009324E2 /)
  double precision, parameter, dimension ( 5 ) :: r4 = (/ &
     0.279195317918525, &
     0.4917317610505968, &
     0.0692910599291889, &
     3.350343815022304, &
     6.012459259764103 /)
  double precision x
  double precision x1
  double precision x2
  double precision y
!
  alngam = 0.0
  ifault = 0
!
!  X <= 0.0
!
  if ( x <= 0.0 ) then

    ifault = 1
    write ( *, * ) ' '
    write ( *, * ) 'ALNGAM - Fatal error!'
    write ( *, * ) '  X <= 0.'
    stop
!
!  0 < X < 0.5
!
  else if ( x < 0.5 ) then

    alngam = - log ( x )
    y = x + 1.0

    if ( y == 1.0 ) then
      return
    end if

    alngam = alngam + x * (((( &
           r1(5)   * y &
         + r1(4) ) * y &
         + r1(3) ) * y &
         + r1(2) ) * y &
         + r1(1) ) / (((( &
                     y &
         + r1(9) ) * y &
         + r1(8) ) * y &
         + r1(7) ) * y &
         + r1(6) )
!
!  0.5 <= X < 1.5
!
  else if ( x < 1.5 ) then

    alngam = ( ( x - 0.5 ) - 0.5 ) * (((( &
           r1(5)   * x &
         + r1(4) ) * x &
         + r1(3) ) * x &
         + r1(2) ) * x &
         + r1(1) ) / (((( &
                     x &
         + r1(9) ) * x &
         + r1(8) ) * x &
         + r1(7) ) * x &
         + r1(6) )
!
!  1.5 <= X < 4.0.
!
  else if ( x < 4.0 ) then

    y = ( x - 1.0 ) - 1.0

    alngam = y * (((( &
           r2(5)   * x &
         + r2(4) ) * x &
         + r2(3) ) * x &
         + r2(2) ) * x &
         + r2(1) ) / (((( &
                     x &
         + r2(9) ) * x &
         + r2(8) ) * x &
         + r2(7) ) * x &
         + r2(6) )
!
!  4.0 <= X < 12.0.
!
  else if ( x < 12.0 ) then

    alngam = (((( &
           r3(5)   * x &
         + r3(4) ) * x &
         + r3(3) ) * x &
         + r3(2) ) * x &
         + r3(1) ) / (((( &
                     x &
         + r3(9) ) * x &
         + r3(8) ) * x &
         + r3(7) ) * x &
         + r3(6) )
!
!  12.0 <= X <= XLGE
!
  else if ( x <= xlge ) then

    y = log ( x )
    alngam = x * ( y - 1.0 ) - 0.5 * y + alr2pi

    x1 = 1.0 / x
    x2 = x1**2

    alngam = alngam + x1 * ( ( &
              r4(3)   * &
         x2 + r4(2) ) * &
         x2 + r4(1) ) / ( ( &
         x2 + r4(5) ) * &
         x2 + r4(4) )
!
!  XLGE < X < XLGST
!
  else if ( x < xlgst ) then

    y = log ( x )
    alngam = x * ( y - 1.0 ) - 0.5 * y + alr2pi
!
!  XLGST <= X
!
  else if ( x >= xlgst ) then

    ifault = 2
    write ( *, * ) ' '
    write ( *, * ) 'ALNGAM - Fatal error!'
    write ( *, * ) '  X is too large.'
    stop

  end if

  return
end
function alnorm ( x, upper )
!
!*******************************************************************************
!
!! ALNORM computes the cumulative density of the standard normal distribution.
!
!
!  Reference:
!
!    I Hill,
!    The Normal Integral,
!    Algorithm AS 66,
!    Applied Statistics,
!    Volume 22, Number 3, pages 424-427, 1973.
!
!  Modified:
!
!    28 March 1999
!
!  Parameters:
!
!    Input, real X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
  double precision, parameter :: a1 = 5.75885480458
  double precision, parameter :: a2 = 2.62433121679
  double precision, parameter :: a3 = 5.92885724438
  double precision, parameter :: b1 = -29.8213557807
  double precision, parameter :: b2 = 48.6959930692
  double precision, parameter :: c1 = -0.000000038052
  double precision, parameter :: c2 = 0.000398064794
  double precision, parameter :: c3 = -0.151679116635
  double precision, parameter :: c4 = 4.8385912808
  double precision, parameter :: c5 = 0.742380924027
  double precision, parameter :: c6 = 3.99019417011
  double precision, parameter :: con = 1.28
  double precision, parameter :: d1 = 1.00000615302
  double precision, parameter :: d2 = 1.98615381364
  double precision, parameter :: d3 = 5.29330324926
  double precision, parameter :: d4 = -15.1508972451
  double precision, parameter :: d5 = 30.789933034
  double precision, parameter :: ltone = 7.0
  double precision, parameter :: p = 0.398942280444
  double precision, parameter :: q = 0.39990348504
  double precision, parameter :: r = 0.398942280385
  double precision, parameter :: utzero = 18.66
!
  double precision alnorm
  logical up
  logical upper
  double precision x
  double precision y
  double precision z
!
  up = upper
  z = x

  if ( z < 0.0 ) then
    up = .not. up
    z = - z
  end if
!
!  TAKE ANOTHER LOOK AT THIS SET OF CONDITIONS.
!
  if ( z > ltone .and. ( ( .not. up ) .or. z > utzero ) ) then

    if ( up ) then
      alnorm = 0.0
    else
      alnorm = 1.0
    end if

    return

  end if

  y = 0.5 * z**2

  if ( z <= con ) then

    alnorm = 0.5 - z * ( p - q * y &
         / ( y + a1 + b1 &
         / ( y + a2 + b2 &
         / ( y + a3 ))))

  else

    alnorm = r * exp ( - y ) &
         / ( z + c1 + d1 &
         / ( z + c2 + d2 &
         / ( z + c3 + d3 &
         / ( z + c4 + d4 &
         / ( z + c5 + d5 &
         / ( z + c6 ))))))

  end if

  if ( .not. up ) then
    alnorm = 1.0 - alnorm
  end if

  return
end
function alogam ( x, ifault )
!
!*******************************************************************************
!
!! ALOGAM computes the logarithm of the Gamma function.
!
!
!  Reference:
!
!    M C Pike and I D Hill,
!    Algorithm 291, Logarithm of Gamma Function,
!    Communications of the ACM,
!    September 1966, volume 9, page 684.
!
!  Modified:
!
!    28 March 1999
!
!  Parameters:
!
!    Input, real X, the argument of the Gamma function.
!    X should be greater than 0.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    1, X <= 0.
!
!    Output, real ALOGAM, the logarithm of the Gamma function of X.
!
  double precision alogam
  double precision f
  integer ifault
  double precision x
  double precision y
  double precision z
!
  if ( x <= 0.0 ) then
    ifault = 1
    alogam = 0.0
    write ( *, * ) ' '
    write ( *, * ) 'ALOGAM - Fatal error!'
    write ( *, * ) '  X <= 0.'
    stop
  end if

  ifault = 0
  y = x

  if ( x < 7.0 ) then

    f = 1.0
    z = y

10      continue

    if ( z < 7.0 ) then
      f = f * z
      z = z + 1.0
      go to 10
    end if

    y = z
    f = - log ( f )

  else

    f = 0.0

  end if

  z = 1.0 / y**2
    
  alogam = f + ( y - 0.5 ) * log ( y ) - y &
       + 0.918938533204673 + &
       ((( &
       - 0.000595238095238   * z &
       + 0.000793650793651 ) * z &
       - 0.002777777777778 ) * z &
       + 0.083333333333333 ) / y

  return
end
function digamma ( x )
!
!*******************************************************************************
!
!! DIGAMMA calculates the digamma or Psi function = d ( LOG ( GAMMA ( X ) ) ) / dX
!
!
!  Reference:
!
!    J Bernardo,
!    Psi ( Digamma ) Function,
!    Algorithm AS 103,
!    Applied Statistics,
!    Volume 25, Number 3, pages 315-317, 1976.
!
!  Modified:
!
!    03 January 2000
!
!  Parameters:
!
!    Input, real X, the argument of the digamma function.
!    0 < X.
!
!    Output, real DIGAMMA, the value of the digamma function at X.
!
  double precision, parameter :: c = 8.5
  double precision, parameter :: d1 = -0.5772156649
  double precision, parameter :: s = 0.00001
  double precision, parameter :: s3 = 0.08333333333
  double precision, parameter :: s4 = 0.0083333333333
  double precision, parameter :: s5 = 0.003968253968
!
  double precision digamma
  double precision r
  double precision x
  double precision y
!
!  The argument must be positive.
!
  if ( x <= 0.0 ) then

    digamma = 0.0
    write ( *, * ) ' '
    write ( *, * ) 'DIGAMMA - Fatal error!'
    write ( *, * ) '  X <= 0.'
    stop
!
!  Use approximation if argument <= S.
!
  else if ( x <= s ) then

    digamma = d1 - 1.0 / x
!
!  Reduce the argument to DIGAMMA(X + N) where (X + N) >= C.
!
  else

    digamma = 0.0
    y = x

10      continue

    if ( y < c ) then
      digamma = digamma - 1.0 / y
      y = y + 1.0
      go to 10
    end if
!
!  Use Stirling's (actually de Moivre's) expansion if argument > C.
!
    r = 1.0 / y**2
    digamma = digamma + log ( y ) - 0.5 / y - r * ( s3 - r * ( s4 - r * s5 ) )

  end if

  return
end
subroutine dirichlet_check ( n, a )
!
!*******************************************************************************
!
!! DIRICHLET_CHECK checks the parameters of the Dirichlet PDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components.
!
!    Input, real A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be positive.
!
  integer n
!
  double precision a(n)
  integer i
  logical positive
!
  positive = .false.

  do i = 1, n

    if ( a(i) < 0.0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'DIRICHLET_CHECK - Fatal error!'
      write ( *, * ) '  A(I) < 0.'
      write ( *, * ) '  For I = ', i
      write ( *, * ) '  A(I) = ', a(i)
      stop
    else if ( a(i) > 0.0 ) then
      positive = .true.
    end if

  end do

  if ( .not. positive ) then
    write ( *, * ) ' '
    write ( *, * ) 'DIRICHLET_CHECK - Fatal error!'
    write ( *, * ) '  All parameters are zero!'
    stop
  end if

  return
end
subroutine dirichlet_estimate ( k, n, x, ix, init, alpha, rlogl, v, g, &
  niter, s, eps, work, ifault )
!
!*******************************************************************************
!
!! DIRICHLET_ESTIMATE estimates the parameters of a Dirichlet distribution.
!
!
!  Reference:
!
!    A. Naryanan,
!    Algorithm AS 266,
!    Maximum Likelihood Estimation of the Parameters of the
!      Dirichlet Distribution,
!    Applied Statistics,
!    Volume 40, Number 2, pages 365-374, 1991.
!
!  Auxilliary routines:
!
!    ALOGAM (CACM algorithm 291 or AS 245),
!    DIGAMA (AS 103), 
!    GAMMAD (AS 239), 
!    PPCHI2 (AS 91), 
!    TRIGAM (AS 121).
!
!  Modified:
!
!    16 June 1999
!
!  Parameters:
!
!    Input, integer K, the number of parameters.
!    2 <= K.
!
!    Input, integer N, the number of observations.
!    K < N.
!
!    Input, real X(IX,K), contains the N by K array of samples
!    from the distribution.  X(I,J) is the J-th component of
!    the I-th sample.
!
!    Input, integer IX, the leading dimension of the array X.
!    N <= IX.
!
!    Input, integer INIT, specifies how the parameter estimates
!    are to be initialized:
!    1, use the method of moments;
!    2, initialize each ALPHA to the minimum of X;
!    otherwise, the input values of ALPHA already contain estimates.
!
!    Input/output, real ALPHA(K).
!    On input, if INIT is not 1 or 2, then ALPHA must contain
!    initial estimates for the parameters.
!    On output, with IFAULT = 0, ALPHA contains the computed
!    estimates for the parameters.
!
!    Output, real RLOGL, the value of the log-likelihood function
!    at the solution point.
!
!    Output, real V(K*(K+1)/2); V(J*(J-1)/2+I) contains the covariance 
!    between ALPHA(I) and ALPHA(J), for I = 1 to J, J = 1 to K.
!
!    Output, real G(K), contains an estimate of the derivative of
!    the log-likelihood with respect to each component of ALPHA.
!
!    Output, integer NITER, contains the number of Newton-Raphson
!    iterations performed.
!
!    Output, real S, the value of the chi-squared statistic.
!
!    Output, real EPS, contains the probability that the chi-squared
!    statistic is less than S.
!
!    Workspace, real WORK(2*K).
!
!    Output, integer IFAULT, error indicator.
!    0, no error, the results were computed successfully;
!    1, K < 2;
!    2, N <= K;
!    3, IX < N;
!    4, if X(I,J) <= 0 for any I or J, or if
!       ABS ( Sum ( 1 <= J <= K ) X(I,J) - 1 ) >= GAMMA = 0.001;
!    5, if IFAULT is returned nonzero from the chi-square 
!       routine PPCHI2;
!    6, if ALPHA(J) <= 0 for any J during any step of the iteration;
!    7, if MAXIT iterations were carried out but convergence
!       was not achieved.
!
  double precision, parameter :: alpha_min = 0.00001
  double precision, parameter :: gamma = 0.0001
  integer, parameter :: it_max = 100
!
  integer ix
  integer k
!
  double precision alngam
  double precision alpha(k)
  double precision an
  double precision beta
  double precision chi2
  double precision digamma
  double precision eps
  double precision g(k)
  double precision gammad
  integer i
  integer i2
  integer ifault
  integer ifault2
  integer init
  integer it_num
  integer j
  integer kk
  integer n
  integer niter
  double precision ppchi2
  double precision rk
  double precision rlogl
  double precision s
  double precision sum
  double precision sum1
  double precision temp
  double precision trigamma
  double precision v((k*(k+1))/2)
  double precision varp1
  double precision work(2*k)
  double precision x(ix,k)
  double precision x_min
  double precision x11
  double precision x12
!
  ifault2 = 0
!
!  Check the input arguments.
!
  if ( k < 2 ) then
    ifault = 1
    write ( *, * ) ' '
    write ( *, * ) 'DIRICHLET_ESTIMATE - Fatal error!'
    write ( *, * ) '  K < 2.'
    stop
  end if

  if ( n <= k ) then
    ifault = 2
    write ( *, * ) ' '
    write ( *, * ) 'DIRICHLET_ESTIMATE - Fatal error!'
    write ( *, * ) '  N <= K.'
    stop
  end if

  if ( ix < n ) then
    ifault = 3
    write ( *, * ) ' '
    write ( *, * ) 'DIRICHLET_ESTIMATE - Fatal error!'
    write ( *, * ) '  IX < N.'
    stop
  end if

  do i = 1, n

    do j = 1, k
      if ( x(i,j) <= 0.0 ) then
        niter = i
        ifault = 4
        write ( *, * ) ' '
        write ( *, * ) 'DIRICHLET_ESTIMATE - Fatal error!'
        write ( *, * ) '  X(I,J) <= 0.'
        stop
      end if
    end do

    sum = 0.0
    do j = 1, k
      sum = sum + x(i,j)
    end do

    if ( abs ( sum - 1.0 ) >= gamma ) then
      ifault = 4
      niter = i
      write ( *, * ) ' '
      write ( *, * ) 'DIRICHLET_ESTIMATE - Fatal error!'
      write ( *, * ) '  Row I does not sum to 1.'
      stop
    end if

  end do

  ifault = 0

  an = dble ( n )
  rk = dble ( k )
  niter = 0
!
!  Calculate initial estimates using the method of moments.
!
  if ( init == 1 ) then

    do j = 1, k - 1
      x12 = 0.0
      do i = 1, n
        x12 = x12 + x(i,j)
      end do
      alpha(j) = x12 / an
    end do

    call rvec_sum ( k-1, alpha, sum )
    alpha(k) = 1.0 - sum

    x12 = 0.0
    do i = 1, n
      x12 = x12 + x(i,1) ** 2
    end do

    x12 = x12 / an
    varp1 = x12 - alpha(1) ** 2

    x11 = ( alpha(1) - x12 ) / varp1
    alpha(1:k) = x11 * alpha(1:k)
!
!  Calculate initial estimates using Ronning's suggestion.
!
  else if ( init == 2 ) then

    x_min = x(1,1)
    do j = 1, k
      do i = 1, n
        x_min = min ( x_min, x(i,j) )
      end do
    end do

    x_min = max ( x_min, alpha_min )

    alpha(1:k) = x_min

  end if
!
!  Check whether any ALPHA's are negative or zero.
!
  do j = 1, k
    if ( alpha(j) <= 0.0 ) then
      ifault = 6
      write ( *, * ) ' '
      write ( *, * ) 'DIRICHLET_ESTIMATE - Fatal error!'
      write ( *, * ) '  For J = ', j
      write ( *, * ) '  ALPHA(J) = ', alpha(j)
      write ( *, * ) '  but ALPHA(J) must be positive.'
      stop
    end if
  end do
!
!  Calculate N * log ( G(J) ) for J = 1,2,...,K and store in WORK array.
!
  do j = 1, k
    sum = 0.0
    do i = 1, n
      sum = sum + log ( x(i,j) )
    end do
    work(j) = sum
  end do
!
!  Call Algorithm AS 91 to compute CHI2, the chi-squared value.
!
  chi2 = ppchi2 ( gamma, rk, ifault )

  if ( ifault > 0 ) then
    ifault = 5
    write ( *, * ) ' '
    write ( *, * ) 'DIRICHLET_ESTIMATE - Fatal error!'
    write ( *, * ) '  PPCHI2 returns error code.'
    stop
  end if
!
!  Carry out the Newton iteration.
!
  do it_num = 1, it_max

    call rvec_sum ( k, alpha, sum )

    sum1 = 0.0
    do j = 1, k
      work(k+j) = trigamma ( alpha(j) )
      sum1 = sum1 + 1.0 / work(k+j)
    end do

    beta = trigamma ( sum )
    beta = an * beta / ( 1.0 - beta * sum1 )

    temp = digamma ( sum )

    do j = 1, k
      g(j) = an * ( temp - digamma ( alpha(j) ) ) + work(j)
    end do
!
!  Calculate the lower triangle of the Variance-Covariance matrix V.
!
    sum = beta / an**2
    do i = 1, k
      do j = 1, i
        kk = j + ( i * ( i - 1 ) ) / 2
        v(kk) = sum / ( work(k+i) * work(k+j) )
        if ( i == j ) then
          v(kk) = v(kk) + 1.0 / ( an * work(k+j) )
        end if
      end do
    end do
!
!  Postmultiply the Variance-Covariance matrix V by G and store
!  in the last K elements of WORK.
!
    do i = 1, k

      sum = 0.0
      i2 = ( i * ( i - 1 ) ) / 2
      do j = 1, i - 1
        sum = sum + v(i2+j) * g(j)
      end do
      do j = i + 1, k
        sum = sum + v(i+(j*(j-1))/2) * g(j)
      end do

      work(k+i) = sum + v((i*(i+1))/2) * g(i)

    end do
!
!  Update the ALPHA's.
!
    niter = it_num

    do j = 1, k
      alpha(j) = alpha(j) + work(k+j)
      alpha(j) = max ( alpha(j), alpha_min )
    end do

    do j = 1, k
      if ( alpha(j) <= 0.0 ) then
        ifault = 6
        write ( *, * ) ' '
        write ( *, * ) 'DIRICHLET_ESTIMATE - Fatal error!'
        write ( *, * ) '  Newton iteration ', it_num
        write ( *, * ) '  Computed ALPHA(J) <= 0.'
        write ( *, * ) '  J = ', j
        write ( *, * ) '  ALPHA(J) = ', alpha(j)
        stop
      end if
    end do
!
!  Test for convergence.
!
    s = 0.0
    do j = 1, k
      s = s + g(j) * work(k+j)
    end do

    if ( s < chi2 ) then
      eps = gammad ( s / 2.0, rk / 2.0, ifault2 )

      call rvec_sum ( k, alpha, sum )

      rlogl = 0.0
      do j = 1, k
        rlogl = rlogl + ( alpha(j) - 1.0 ) * work(j) - an * &
          alngam ( alpha(j), ifault2 )
      end do

      rlogl = rlogl + an * alngam ( sum, ifault2 )

      return

    end if

  end do

  ifault = 7

  write ( *, * ) ' '
  write ( *, * ) 'DIRICHLET_ESTIMATE - Fatal error!'
  write ( *, * ) '  No convergence.'

  stop
end
subroutine dirichlet_mean ( n, a, mean )
!
!*******************************************************************************
!
!! DIRICHLET_MEAN returns the means of the Dirichlet PDF.
!
!
!  Modified:
!
!    23 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components.
!
!    Input, real A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be positive.
!
!    Output, real MEAN(N), the means of the PDF.
!
  integer n
!
  double precision a(n)
  double precision a_sum
  integer i
  double precision mean(n)
!
  call dirichlet_check ( n, a )

  call rvec_sum ( n, a, a_sum )

  mean(1:n) = a(1:n) / a_sum

  return
end
subroutine dirichlet_mix_mean ( comp_max, comp_num, elem_num, a, comp_weight, &
  mean )
!
!*******************************************************************************
!
!! DIRICHLET_MIX_MEAN returns the means of a Dirichlet mixture PDF.
!
!
!  Modified:
!
!    29 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer COMP_MAX, the leading dimension of A, which must
!    be at least COMP_NUM.
!
!    Input, integer COMP_NUM, the number of components in the Dirichlet
!    mixture density, that is, the number of distinct Dirichlet PDF's
!    that are mixed together.
!
!    Input, integer ELEM_NUM, the number of elements of an observation.
!
!    Input, real A(COMP_MAX,ELEM_NUM), the probabilities for element ELEM_NUM
!    in component COMP_NUM.
!    Each A(I,J) should be greater than or equal to 0.0.
!
!    Input, integer COMP_WEIGHT(COMP_NUM), the mixture weights of the densities.
!    These do not need to be normalized.  The weight of a given component is
!    the relative probability that that component will be used to generate
!    the sample.
!
!    Output, real MEAN(ELEM_NUM), the means for each element.
!
  integer comp_max
  integer comp_num
  integer elem_num
!
  double precision a(comp_max,elem_num)
  double precision a_sum(comp_num)
  integer comp_i
  double precision comp_weight(comp_num)
  double precision comp_weight_sum
  integer elem_i
  double precision mean(elem_num)
!
!  Check.
!
  do comp_i = 1, comp_num

    do elem_i = 1, elem_num
      if ( a(comp_i,elem_i) < 0.0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'DIRICHLET_MIX_MEAN - Fatal error!'
        write ( *, * ) '  A(COMP,ELEM) < 0.'
        write ( *, * ) '  COMP = ', comp_i
        write ( *, * ) '  ELEM = ', elem_i
        write ( *, * ) '  A(COMP,ELEM) = ', a(comp_i,elem_i)
        stop
      end if
    end do

  end do

  call rvec_sum ( comp_num, comp_weight, comp_weight_sum )

  do comp_i = 1, comp_num
    a_sum(comp_i) = 0.0
    do elem_i = 1, elem_num
      a_sum(comp_i) = a_sum(comp_i) + a(comp_i,elem_i)
    end do
  end do

  do elem_i = 1, elem_num
    mean(elem_i) = 0.0
    do comp_i = 1, comp_num
      mean(elem_i) = mean(elem_i) + comp_weight(comp_i) * a(comp_i,elem_i) &
        / a_sum(comp_i)
    end do
    mean(elem_i) = mean(elem_i) / comp_weight_sum
  end do

  return
end
subroutine dirichlet_mix_sample ( comp_max, comp_num, elem_num, a, &
  comp_weight, iseed, comp, x )
!
!*******************************************************************************
!
!! DIRICHLET_MIX_SAMPLE samples a Dirichlet mixture PDF.
!
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer COMP_MAX, the leading dimension of A, which must
!    be at least COMP_NUM.
!
!    Input, integer COMP_NUM, the number of components in the Dirichlet 
!    mixture density, that is, the number of distinct Dirichlet PDF's
!    that are mixed together.
!
!    Input, integer ELEM_NUM, the number of elements of an observation.
!
!    Input, real A(COMP_MAX,ELEM_NUM), the probabilities for element ELEM_NUM
!    in component COMP_NUM.  
!    Each A(I,J) should be greater than or equal to 0.0.
!
!    Input, integer COMP_WEIGHT(COMP_NUM), the mixture weights of the densities.
!    These do not need to be normalized.  The weight of a given component is
!    the relative probability that that component will be used to generate
!    the sample.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, integer COMP, the index of the component of the Dirichlet
!    mixture that was chosen to generate the sample.
!
!    Output, real X(ELEM_NUM), a sample of the PDF.
!
  integer comp_max
  integer comp_num
  integer elem_num
!
  double precision a(comp_max,elem_num)
  integer comp
  integer comp_i
  double precision comp_weight(comp_num)
  double precision comp_weight_sum
  integer elem_i
  integer iseed
  double precision r
  double precision sum
  double precision x(elem_num)
  double precision x_sum
!
!  Check.
!
  do comp_i = 1, comp_num

    do elem_i = 1, elem_num
      if ( a(comp_i,elem_i) < 0.0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'DIRICHLET_MIX_SAMPLE - Fatal error!'
        write ( *, * ) '  A(COMP,ELEM) < 0.'
        write ( *, * ) '  COMP = ', comp_i
        write ( *, * ) '  ELEM = ', elem_i
        write ( *, * ) '  A(COMP,ELEM) = ', a(comp_i,elem_i)
        stop
      end if
    end do

  end do
!
!  Choose a particular density MIX.
!
  call rvec_sum ( comp_num, comp_weight, comp_weight_sum )

  call r_random ( r, 0.d0, comp_weight_sum, iseed )

  comp = 1
  sum = 0.0

10    continue

  if ( comp < comp_num ) then

    sum = sum + comp_weight(comp)

    if ( r > sum ) then
      comp = comp + 1
      go to 10
    end if

  end if
!
!  Sample density COMP.
!
  do elem_i = 1, elem_num
    call gamma_sample ( a(comp,elem_i), 1d0, iseed, x(elem_i) )
  end do
!
!  Normalize the result.
!
  call rvec_sum ( elem_num, x, x_sum )

  do elem_i = 1, elem_num
    x(elem_i) = x(elem_i) / x_sum
  end do

  return
end
subroutine dirichlet_sample ( n, a, iseed, x )
!
!*******************************************************************************
!
!! DIRICHLET_SAMPLE samples the Dirichlet PDF.
!
!
!  Reference:
!
!    Jerry Banks, editor,
!    Handbook of Simulation,
!    Engineering and Management Press Books, 1998, page 169.
!
!  Modified:
!
!    23 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components.
!
!    Input, real A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be
!    positive.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X(N), a sample of the PDF.  The entries of X should
!    sum to 1.
!
  integer n
!
  double precision a(n)
  integer i
  integer iseed
  double precision x(n)
  double precision x_sum
!
  call dirichlet_check ( n, a )

  do i = 1, n
    call gamma_sample ( a(i), 1d0, iseed, x(i) )
  end do
!
!  Normalize the result.
!
  call rvec_sum ( n, x, x_sum )

  do i = 1, n
    x(i) = x(i) / x_sum
  end do

  return
end
subroutine dirichlet_variance ( n, a, variance )
!
!*******************************************************************************
!
!! DIRICHLET_VARIANCE returns the variances of the Dirichlet PDF.
!
!
!  Modified:
!
!    03 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components.
!
!    Input, real A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be positive.
!
!    Output, real VARIANCE(N), the variances of the PDF.
!
  integer n
!
  double precision a(n)
  double precision a_sum
  integer i
  double precision variance(n)
!
  call dirichlet_check ( n, a )

  call rvec_sum ( n, a, a_sum )

  do i = 1, n
    variance(i) = a(i) * ( a_sum - a(i) ) / ( a_sum**2 * ( a_sum + 1.0 ) )
  end do

  return
end
subroutine exponential_1_sample ( iseed, x )
!
!*******************************************************************************
!
!! EXPONENTIAL_1_SAMPLE samples the Exponential PDF with parameter 1.
!
!
!  Modified:
!
!    16 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  double precision a
  double precision b
  double precision cdf
  integer iseed
  double precision uniform_01_sample
  double precision x
!
  cdf = uniform_01_sample ( iseed )

  a = 0.0
  b = 1.0
  call exponential_cdf_inv ( cdf, a, b, x )

  return
end
subroutine exponential_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! EXPONENTIAL_CDF_INV inverts the Exponential CDF.
!
!
!  Modified:
!
!    01 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real X, the corresponding argument.
!
  double precision a
  double precision b
  double precision cdf
  double precision x
!
!  Check.
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'EXPONENTIAL_CDF_INV - Fatal error!'
    write ( *, * ) '  B <= 0.0'
    stop
  end if

  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'EXPONENTIAL_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if
!
  x = a - b * log ( 1.0 - cdf )

  return
end
function gamain ( x, p, ifault )
!
!*******************************************************************************
!
!! GAMAIN computes the incomplete Gamma ratio.
!
!
!  Discussion:
!
!    A series expansion is used if P > X or X <= 1.  Otherwise, a
!    continued fraction approximation is used.
!
!  Reference:
!
!    G Bhattacharjee,
!    The Incomplete Gamma Integral,
!    Algorithm AS 32,
!    Applied Statistics,
!    Volume 19, Number 3, pages 285-287, 1970.
!
!  Modified:
!
!    30 March 1999
!
!  Parameters:
!
!    Input, real X, P, the parameters of the incomplete gamma ratio.
!    0 <= X, and 0 < P.
!
!    Output, integer IFAULT, error flag.
!    0, no errors.
!    1, P <= 0.
!    2, X < 0.
!    3, underflow.
!    4, error return from the Log Gamma function.
!
  double precision, parameter :: acu = 1.0E-08
  double precision, parameter :: oflo = 1.0E+37
  double precision, parameter :: uflo = 1.0E-37
!
  double precision a
  double precision alngam
  double precision an
  double precision arg
  double precision b
  double precision dif
  double precision factor
  double precision g
  double precision gamain
  double precision gin
  integer i
  integer ifault
  double precision p
  double precision pn(6)
  double precision rn
  double precision term
  double precision x
!
!  Check the input.
!
  if ( p <= 0.0 ) then
    ifault = 1
    gamain = 0.0
    write ( *, * ) ' '
    write ( *, * ) 'GAMAIN - Fatal error!'
    write ( *, * ) '  P <= 0.'
    stop
  end if

  if ( x < 0.0 ) then
    ifault = 2
    gamain = 0.0
    write ( *, * ) ' '
    write ( *, * ) 'GAMAIN - Fatal error!'
    write ( *, * ) '  X < 0.'
    stop
  end if

  if ( x == 0.0 ) then
    ifault = 0
    gamain = 0.0
    return
  end if

  g = alngam ( p, ifault )

  if ( ifault /= 0 ) then
    ifault = 4
    gamain = 0.0
    write ( *, * ) ' '
    write ( *, * ) 'GAMAIN - Fatal error!'
    write ( *, * ) '  ALNGAM returns error code.'
    stop
  end if

  arg = p * log ( x ) - x - g

  if ( arg < log ( uflo ) ) then
    ifault = 3
    gamain = 0.0
    write ( *, * ) ' '
    write ( *, * ) 'GAMAIN - Fatal error!'
    write ( *, * ) '  Underflow.'
    stop
  end if

  ifault = 0
  factor = exp ( arg )
!
!  Calculation by series expansion.
!
  if ( x <= 1.0 .or. x < p ) then

    gin = 1.0
    term = 1.0
    rn = p

20      continue

    rn = rn + 1.0
    term = term * x / rn
    gin = gin + term

    if ( term > acu ) then
      go to 20
    end if

    gamain = gin * factor / p
    return

  end if
!
!  Calculation by continued fraction.
!
  a = 1.0 - p
  b = a + x + 1.0
  term = 0.0

  pn(1) = 1.0
  pn(2) = x
  pn(3) = x + 1.0
  pn(4) = x * b

  gin = pn(3) / pn(4)

32    continue

  a = a + 1.0
  b = b + 2.0
  term = term + 1.0
  an = a * term
  do i = 1, 2
    pn(i+4) = b * pn(i+2) - an * pn(i)
  end do

  if ( pn(6) /= 0.0 ) then

    rn = pn(5) / pn(6)
    dif = abs ( gin - rn )
!
!  Absolute error tolerance satisfied?
!
    if ( dif <= acu ) then
!
!  Relative error tolerance satisfied?
!
      if ( dif <= acu * rn ) then
        gamain = 1.0 - factor * gin
        return
      end if

    end if

    gin = rn

  end if

  do i = 1, 4
    pn(i) = pn(i+2)
  end do

  if ( abs ( pn(5) ) < oflo ) then
    go to 32
  end if

  do i = 1, 4
    pn(i) = pn(i) / oflo
  end do

  go to 32
end
function gamma_log ( x )
!
!*******************************************************************************
!
!! GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
!
!
!  Discussion:
!
!    Computation is based on an algorithm outlined in references 1 and 2.
!    The program uses rational functions that theoretically approximate
!    LOG(GAMMA(X)) to at least 18 significant decimal digits.  The
!    approximation for X > 12 is from reference 3, while approximations
!    for X < 12.0 are similar to those in reference 1, but are unpublished.
!    The accuracy achieved depends on the arithmetic system, the compiler,
!    intrinsic functions, and proper selection of the machine-dependent
!    constants.
!
!  Modified:
!
!    16 June 1999
!
!  Authors:
!
!    W. J. Cody and L. Stoltz
!    Argonne National Laboratory
!
!  References:
!
!    # 1)
!    W. J. Cody and K. E. Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
!    Math. Comp.
!    Volume 21, 1967, pages 198-203.
!
!    # 2)
!    K. E. Hillstrom,
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
!    May 1969.
!
!    # 3)
!    Hart, Et. Al.,
!    Computer Approximations,
!    Wiley and sons, New York, 1968.
!
!  Parameters:
!
!    Input, real X, the argument of the Gamma function.  X must be positive.
!
!    Output, real GAMMA_LOG, the logarithm of the Gamma function of X.
!    If X <= 0.0, or if overflow would occur, the program returns the
!    value XINF, the largest representable floating point number.
!
!*******************************************************************************
!
!  Explanation of machine-dependent constants
!
!  BETA   - radix for the floating-point representation.
!
!  MAXEXP - the smallest positive power of BETA that overflows.
!
!  XBIG   - largest argument for which LN(GAMMA(X)) is representable
!           in the machine, i.e., the solution to the equation
!             LN(GAMMA(XBIG)) = BETA**MAXEXP.
!
!  XINF   - largest machine representable floating-point number;
!           approximately BETA**MAXEXP.
!
!  EPS    - The smallest positive floating-point number such that
!           1.0+EPS .GT. 1.0
!
!  FRTBIG - Rough estimate of the fourth root of XBIG
!
!
!  Approximate values for some important machines are:
!
!                            BETA      MAXEXP         XBIG
!
!  CRAY-1        (S.P.)        2        8191       9.62E+2461
!  Cyber 180/855
!    under NOS   (S.P.)        2        1070       1.72E+319
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)        2         128       4.08E+36
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)        2        1024       2.55D+305
!  IBM 3033      (D.P.)       16          63       4.29D+73
!  VAX D-Format  (D.P.)        2         127       2.05D+36
!  VAX G-Format  (D.P.)        2        1023       1.28D+305
!
!
!                            XINF        EPS        FRTBIG
!
!  CRAY-1        (S.P.)   5.45E+2465   7.11E-15    3.13E+615
!  Cyber 180/855
!    under NOS   (S.P.)   1.26E+322    3.55E-15    6.44E+79
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)   3.40E+38     1.19E-7     1.42E+9
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)   1.79D+308    2.22D-16    2.25D+76
!  IBM 3033      (D.P.)   7.23D+75     2.22D-16    2.56D+18
!  VAX D-Format  (D.P.)   1.70D+38     1.39D-17    1.20D+9
!  VAX G-Format  (D.P.)   8.98D+307    1.11D-16    1.89D+76
!
  double precision, parameter :: d1 = - 5.772156649015328605195174e-1
  double precision, parameter :: d2 =   4.227843350984671393993777e-1
  double precision, parameter :: d4 =   1.791759469228055000094023e0
  double precision, parameter :: EPS = 1.19e-7
  double precision, parameter :: FRTBIG = 1.42e9
  double precision, parameter :: PNT68 = 0.6796875
  double precision, parameter :: SQRTPI = 0.9189385332046727417803297
  double precision, parameter :: XBIG = 4.08e36
  double precision, parameter :: XINF = 3.401e38
!
  double precision, parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728e-03, &
     8.4171387781295e-04, &
    -5.952379913043012e-04, &
     7.93650793500350248e-04, &
    -2.777777777777681622553e-03, &
     8.333333333333333331554247e-02, &
     5.7083835261e-03 /)
  double precision corr
  integer i
  double precision gamma_log
  double precision, parameter, dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888e0, &
    2.018112620856775083915565e2, &
    2.290838373831346393026739e3, &
    1.131967205903380828685045e4, &
    2.855724635671635335736389e4, &
    3.848496228443793359990269e4, &
    2.637748787624195437963534e4, &
    7.225813979700288197698961e3 /)
  double precision, parameter, dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064e0, &
    5.424138599891070494101986e2, &
    1.550693864978364947665077e4, &
    1.847932904445632425417223e5, &
    1.088204769468828767498470e6, &
    3.338152967987029735917223e6, &
    5.106661678927352456275255e6, &
    3.074109054850539556250927e6 /)
  double precision, parameter, dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062e4, &
    2.426813369486704502836312e6, &
    1.214755574045093227939592e8, &
    2.663432449630976949898078e9, &
    2.940378956634553899906876e10, &
    1.702665737765398868392998e11, &
    4.926125793377430887588120e11, &
    5.606251856223951465078242e11 /)
  double precision, parameter, dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036e1, &
    1.113332393857199323513008e3, &
    7.738757056935398733233834e3, &
    2.763987074403340708898585e4, &
    5.499310206226157329794414e4, &
    6.161122180066002127833352e4, &
    3.635127591501940507276287e4, &
    8.785536302431013170870835e3 /)
  double precision, parameter, dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942e2, &
    7.765049321445005871323047e3, &
    1.331903827966074194402448e5, &
    1.136705821321969608938755e6, &
    5.267964117437946917577538e6, &
    1.346701454311101692290052e7, &
    1.782736530353274213975932e7, &
    9.533095591844353613395747e6 /)
  double precision, parameter, dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843e3, &
    6.393885654300092398984238e5, &
    4.135599930241388052042842e7, &
    1.120872109616147941376570e9, &
    1.488613728678813811542398e10, &
    1.016803586272438228077304e11, &
    3.417476345507377132798597e11, &
    4.463158187419713286462081e11 /)
  double precision res
  double precision x
  double precision xden
  double precision xm1
  double precision xm2
  double precision xm4
  double precision xnum
  double precision xsq
!
!  Return immediately if the argument is out of range.
!
  if ( x <= 0.0 .or. x > XBIG ) then
    gamma_log = XINF
    return
  end if

  if ( x <= EPS ) then

    res = - log ( x )

  else if ( x <= 1.5 ) then

    if ( x < PNT68 ) then
      corr = - log ( x )
      xm1 = x
    else
      corr = 0.0
      xm1 = ( x - 0.5 ) - 0.5
    end if

    if ( x <= 0.5 .or. x >= PNT68 ) then

      xden = 1.0
      xnum = 0.0

      do i = 1, 8
        xnum = xnum * xm1 + p1(i)
        xden = xden * xm1 + q1(i)
      end do

      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

    else

      xm2 = ( x - 0.5 ) - 0.5
      xden = 1.0
      xnum = 0.0
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

    end if

  else if ( x <= 4.0 ) then

    xm2 = x - 2.0
    xden = 1.0
    xnum = 0.0
    do i = 1, 8
      xnum = xnum * xm2 + p2(i)
      xden = xden * xm2 + q2(i)
    end do

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) )

  else if ( x <= 12.0 ) then

    xm4 = x - 4.0
    xden = - 1.0
    xnum = 0.0
    do i = 1, 8
      xnum = xnum * xm4 + p4(i)
      xden = xden * xm4 + q4(i)
    end do

    res = d4 + xm4 * ( xnum / xden )

  else

    res = 0.0

    if ( x <= FRTBIG ) then

      res = c(7)
      xsq = x * x

      do i = 1, 6
        res = res / xsq + c(i)
      end do

    end if

    res = res / x
    corr = log ( x )
    res = res + SQRTPI - 0.5 * corr
    res = res + x * ( corr - 1.0 )

  end if

  gamma_log = res

  return
end
subroutine gamma_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! GAMMA_SAMPLE samples the Gamma PDF.
!
!
!  References:
!
!    J H Ahrens and U Dieter,
!    Generating Gamma Variates by a Modified Rejection Technique,
!    Communications of the ACM,
!    Volume 25, Number 1, January 1982, pages 47 - 54.
!
!    J H Ahrens and U Dieter,
!    Computer Methods for Sampling from Gamma, Beta, Poisson and
!      Binomial Distributions.
!    Computing, Volume 12, 1974, pages 223 - 246.
!
!  Modified:
!
!    13 February 1999
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  double precision, parameter :: a1 =  0.3333333
  double precision, parameter :: a2 = -0.2500030
  double precision, parameter :: a3 =  0.2000062
  double precision, parameter :: a4 = -0.1662921
  double precision, parameter :: a5 =  0.1423657
  double precision, parameter :: a6 = -0.1367177
  double precision, parameter :: a7 =  0.1233795
  double precision, parameter :: e1 = 1.0
  double precision, parameter :: e2 = 0.4999897
  double precision, parameter :: e3 = 0.1668290
  double precision, parameter :: e4 = 0.0407753
  double precision, parameter :: e5 = 0.0102930
  double precision, parameter :: q1 =  0.04166669
  double precision, parameter :: q2 =  0.02083148
  double precision, parameter :: q3 =  0.00801191
  double precision, parameter :: q4 =  0.00144121
  double precision, parameter :: q5 = -0.00007388
  double precision, parameter :: q6 =  0.00024511
  double precision, parameter :: q7 =  0.00024240
!
  double precision a
  double precision b
  double precision bcoef
  double precision c
  double precision d
  double precision e
  integer iseed
  double precision p
  double precision q
  double precision q0
  double precision r
  double precision s
  double precision si
  double precision s2
  double precision t
  double precision u
  double precision uniform_01_sample
  double precision v
  double precision w
  double precision x
!
!  A < 1.
!
  if ( a < 1.0 ) then

10      continue

    p = ( 1.0 + 0.3678794 * a ) * uniform_01_sample ( iseed )

    call exponential_1_sample ( iseed, e )

    if ( p >= 1.0 ) then

      x = - log ( ( 1.0 + 0.3678794 * a - p ) / a )

      if ( e >= ( 1.0 - a ) * log ( x ) ) then
        x = x / b
        return
      end if

    else

      x = exp ( log ( p ) / a )

      if ( e >= x ) then
        x = x / b
        return
      end if

    end if

    go to 10
!
!  1 <= A.
!
  else

    s2 = a - 0.5
    s = sqrt ( a - 0.5 )
    d = sqrt ( 32.0 ) - 12.0 * sqrt ( a - 0.5 )

    call normal_01_sample ( iseed, t )
    x = ( sqrt ( a - 0.5 ) + 0.5 * t )**2

    if ( t >= 0.0 ) then
      x = x / b
      return
    end if

    u = uniform_01_sample ( iseed )

    if ( d * u <= t**3 ) then
      x = x / b
      return
    end if

    r = 1.0 / a
    q0 = ( ( ( ( ( ( &
           q7   * r &
         + q6 ) * r &
         + q5 ) * r &
         + q4 ) * r &
         + q3 ) * r &
         + q2 ) * r &
         + q1 ) * r 

    if ( a <= 3.686 ) then
      bcoef = 0.463 + s - 0.178 * s2
      si = 1.235
      c = 0.195 / s - 0.079 + 0.016 * s
    else if ( a <= 13.022 ) then
      bcoef = 1.654 + 0.0076 * s2
      si = 1.68 / s + 0.275
      c = 0.062 / s + 0.024
    else
      bcoef = 1.77
      si = 0.75
      c = 0.1515 / s
    end if

    if ( sqrt ( a - 0.5 ) + 0.5 * t > 0.0 ) then

      v = 0.5 * t / s

      if ( abs ( v ) > 0.25 ) then
        q = q0 - s * t + 0.25 * t**2 + 2.0 * s2 * log ( 1.0 + v )
      else
        q = q0 + 0.5 * t**2 * ( ( ( ( ( ( &
               a7   * v &
             + a6 ) * v &
             + a5 ) * v &
             + a4 ) * v &
             + a3 ) * v &
             + a2 ) * v &
             + a1 ) * v
      end if

      if ( log ( 1.0 - u ) <= q ) then
        x = x / b
        return
      end if

    end if

20      continue

    call exponential_1_sample ( iseed, e )

    u = 2.0 * uniform_01_sample ( iseed ) - 1.0
    t = bcoef + sign ( si * e, u )

    if ( t >= - 0.7187449 ) then

      v = 0.5 * t / s

      if ( abs ( v ) > 0.25 ) then
        q = q0 - s * t + 0.25 * t**2 + 2.0 * s2 * log ( 1.0 + v )
      else
        q = q0 + 0.5 * t**2 * ( ( ( ( ( ( &
               a7   * v &
             + a6 ) * v &
             + a5 ) * v &
             + a4 ) * v &
             + a3 ) * v &
             + a2 ) * v &
             + a1 ) * v
      end if

      if ( q > 0.0 ) then

        if ( q > 0.5 ) then
          w = exp ( q ) - 1.0
        else
          w = ( ( ( ( &
                 e5   * q &
               + e4 ) * q &
               + e3 ) * q &
               + e2 ) * q &
               + e1 ) * q
        end if

        if ( c * abs ( u ) <= w * exp ( e - 0.5 * t**2 ) ) then
          x = ( s + 0.5 * t )**2 / b
          return
        end if

      end if

    end if

    go to 20

  end if

end
function gammad ( x, p, ifault )
!
!*******************************************************************************
!
!! GAMMAD computes the Incomplete Gamma Integral
!
!
!  Reference:
!
!    B Shea,
!    Chi-squared and Incomplete Gamma Integral,
!    Algorithm AS 239,
!    Applied Statistics,
!    Volume 37, Number 3, pages 466-473, 1988.
!
!  Auxiliary functions:
!
!    ALOGAM = logarithm of the gamma function, 
!    ALNORM = algorithm AS66
!
!  Modified:
!
!    31 March 1999
!
!  Parameters:
!
!    Input, real X, P, the parameters of the incomplete gamma ratio.
!    0 <= X, and 0 < P.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    1, X < 0 or P <= 0.
!
!    Output, real GAMMAD, the value of the incomplete Gamma integral.
!
  double precision, parameter :: elimit = - 88.0
  double precision, parameter :: oflo = 1.0E+37
  double precision, parameter :: plimit = 1000.0
  double precision, parameter :: tol = 1.0E-14
  double precision, parameter :: xbig = 1.0E+8
!
  double precision a
  double precision alnorm
  double precision alngam
  double precision an
  double precision arg
  double precision b
  double precision c
  double precision gammad
  integer ifault
  double precision p
  double precision pn1
  double precision pn2
  double precision pn3
  double precision pn4
  double precision pn5
  double precision pn6
  double precision rn
  logical upper
  double precision x
!
  gammad = 0.0
!
!  Check the input.
!
  if ( x < 0.0 ) then
    ifault = 1
    write ( *, * ) ' '
    write ( *, * ) 'GAMMAD - Fatal error!'
    write ( *, * ) '  X < 0.'
    stop
  end if

  if ( p <= 0.0 ) then
    ifault = 1
    write ( *, * ) ' '
    write ( *, * ) 'GAMMAD - Fatal error!'
    write ( *, * ) '  P <= 0.'
    stop
  end if

  ifault = 0

  if ( x == 0.0 ) then
    gammad = 0.0
    return
  end if
!
!  If P is large, use a normal approximation.
!
  if ( p > plimit ) then

    pn1 = 3.0 * sqrt ( p ) * ( ( x / p )**( 1.0 / 3.0 ) + 1.0 / ( 9.0 * p ) &
      - 1.0 )

    upper = .false.
    gammad = alnorm ( pn1, upper )
    return

  end if
!
!  If X is large set GAMMAD = 1.
!
  if ( x > xbig ) then
    gammad = 1.0
    return
  end if
!
!  Use Pearson's series expansion.
!  (Note that P is not large enough to force overflow in ALOGAM).
!  No need to test IFAULT on exit since P > 0.
!
  if ( x <= 1.0 .or. x < p ) then

    arg = p * log ( x ) - x - alngam ( p + 1.0, ifault )
    c = 1.0
    gammad = 1.0
    a = p

   40   continue

    a = a + 1.0
    c = c * x / a
    gammad = gammad + c

    if ( c > tol ) then
      go to 40
    end if

    arg = arg + log ( gammad )

    if ( arg >= elimit ) then
      gammad = exp ( arg )
    else
      gammad = 0.0
    end if
!
!  Use a continued fraction expansion.
!
  else 

    arg = p * log ( x ) - x - alngam ( p, ifault )
    a = 1.0 - p
    b = a + x + 1.0
    c = 0.0
    pn1 = 1.0
    pn2 = x
    pn3 = x + 1.0
    pn4 = x * b
    gammad = pn3 / pn4

   60   continue

    a = a + 1.0
    b = b + 2.0
    c = c + 1.0
    an = a * c
    pn5 = b * pn3 - an * pn1
    pn6 = b * pn4 - an * pn2

    if ( pn6 /= 0.0 ) then

      rn = pn5 / pn6

      if ( abs ( gammad - rn ) <= min ( tol, tol * rn ) ) then
        go to 80
      end if

      gammad = rn

    end if

    pn1 = pn3
    pn2 = pn4
    pn3 = pn5
    pn4 = pn6
!
!  Re-scale terms in continued fraction if terms are large.
!
    if ( abs ( pn5 ) >= oflo ) then
      pn1 = pn1 / oflo
      pn2 = pn2 / oflo
      pn3 = pn3 / oflo
      pn4 = pn4 / oflo
    end if

    go to 60

   80   continue

    arg = arg + log ( gammad )

    if ( arg >= elimit ) then
      gammad = 1.0 - exp ( arg )
    else
      gammad = 1.0
    end if

  end if

  return
end
function gammds ( x, p, ifault )
!
!*******************************************************************************
!
!! GAMMDS computes the incomplete Gamma integral.
!
!
!  Discussion:
!
!    The parameters must be positive.  An infinite series is used.
!
!  Reference:
!
!    Chi_Leung Lau,
!    A Simple Series for the Incomplete Gamma Integral,
!    Algorithm AS 147,
!    Applied Statistics,
!    Volume 29, Number 1, pages 113-114, 1980.
!
!  Auxiliary function:
!
!    ALNGAM = CACM algorithm 291
!
!  Modified:
!
!    30 March 1999
!
!  Parameters:
!
!    Input, real X, P, the arguments of the incomplete Gamma integral.
!    X and P must be greater than 0.
!
!    Output, integer IFAULT, error flag.
!    0, no errors.
!    1, X <= 0 or P <= 0.
!    2, underflow during the computation.
!
!    Output, real GAMMDS, the value of the incomplete Gamma integral.
!
  double precision, parameter :: e = 1.0E-09
  double precision, parameter :: uflo = 1.0E-37
!
  double precision a
  double precision alngam
  double precision arg
  double precision c
  double precision f
  double precision gammds
  integer ifault
  integer ifault2
  double precision p
  double precision x
!
!  Check the input.
!
  if ( x <= 0.0 ) then
    ifault = 1
    gammds = 0.0
    write ( *, * ) ' '
    write ( *, * ) 'GAMMDS - Fatal error!'
    write ( *, * ) '  X <= 0.'
    stop
  end if

  if ( p <= 0.0 ) then
    ifault = 1
    gammds = 0.0
    write ( *, * ) ' '
    write ( *, * ) 'GAMMDS - Fatal error!'
    write ( *, * ) '  P <= 0.'
    stop
  end if
!
!  ALNGAM is the natural logarithm of the gamma function.
!
  ifault2 = 0
  arg = p * log ( x ) - alngam ( p + 1.0, ifault2 ) - x

  if ( arg < log ( uflo ) ) then
    gammds = 0.0
    ifault = 2
    write ( *, * ) ' '
    write ( *, * ) 'GAMMDS - Warning!'
    write ( *, * ) '  Underflow.'
    return
  end if

  f = exp ( arg )

  if ( f == 0.0 ) then
    gammds = 0.0
    ifault = 2
    write ( *, * ) ' '
    write ( *, * ) 'GAMMAD - Warning!'
    write ( *, * ) '  Underflow.'
    return
  end if

  ifault = 0
!
!  Series begins.
!
  c = 1.0
  gammds = 1.0
  a = p

10    continue

  a = a + 1.0
  c = c * x / a
  gammds = gammds + c

  if ( c > e * gammds ) then
    go to 10
  end if

  gammds = gammds * f

  return
end
subroutine get_seed ( iseed )
!
!*******************************************************************************
!
!! GET_SEED returns a seed for the random number generator.
!
!
!  Discussion:
!
!    The seed depends on the current time, and ought to be (slightly)
!    different every millisecond.  Once the seed is obtained, a random
!    number generator should be called a few times to further process
!    the seed.
!
!    If a FORTRAN77 compiler is used, then a different set of routines
!    must be invoked in order to get the date and time.
!
!  Modified:
!
!    23 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ISEED, a pseudorandom seed value.
!
  integer, parameter :: I_MAX = 2147483647
!
  integer iseed
  double precision temp
  integer values(8)
!
!  FORTRAN90 declarations
!
  character ( len = 10 ) time90
  character ( len = 8 ) today90
  character ( len = 5 ) zone
!
!  FORTRAN77 declarations
!
!     integer day
!     integer hour
!     integer minute
!     integer month
!     integer second
!     character*8 time77
!     character*9 today77
!     integer year
!
!  FORTRAN90
!
  call date_and_time ( today90, time90, zone, values )
!
!  FORTRAN77
!
!     call date ( today77 )
!     call s_to_ymd ( today77, 'DD-NNN-YY', year, month, day )
!     call time ( time77 )
!     call s_to_hms ( time77, 'hh:mm:ss', hour, minute, second )
!
!     if ( year == 99 ) then
!       year = 1900 + year
!     else
!       year = 2000 + year
!     end if
!
!     values(1) = year
!     values(2) = month
!     values(3) = day
!     values(4) = 0
!     values(5) = hour
!     values(6) = minute
!     values(7) = second
!     values(8) = 0
!
  temp = 0.0

  temp = temp + dble ( values(2) - 1 ) / 11.0
  temp = temp + dble ( values(3) - 1 ) / 30.0
  temp = temp + dble ( values(5) ) / 23.0
  temp = temp + dble ( values(6) ) / 59.0
  temp = temp + dble ( values(7) ) / 59.0
  temp = temp + dble ( values(8) ) / 999.0
  temp = temp / 6.0

  if ( temp <= 0.0 ) then
    temp = 1.0 / 3.0
  else if ( temp >= 1.0 ) then
    temp = 2.0 / 3.0
  end if

  iseed = int ( dble ( I_MAX ) * temp )
!
!  Never use a seed of 0 or I_MAX.
!
  if ( iseed == 0 ) then
    iseed = 1
  end if

  if ( iseed == I_MAX ) then
    iseed = I_MAX - 1
  end if

  return
end
function lngamma ( z, ier )
!
!*******************************************************************************
!
!! LNGAMMA computes Log(Gamma(X)) using a Lanczos approximation.
!
!
!  Discussion:
!
!    This algorithm is not part of the Applied Statistics algorithms.   
!    It is slower but gives 14 or more significant decimal digits 
!    accuracy, except around X = 1 and X = 2.   The Lanczos series from 
!    which this algorithm is derived is interesting in that it is a 
!    convergent series approximation for the gamma function, whereas 
!    the familiar series due to De Moivre (and usually wrongly called 
!    Stirling's approximation) is only an asymptotic approximation, as
!    is the true and preferable approximation due to Stirling.
!
!  Reference:
!
!    C. Lanczos,
!    A precision approximation of the gamma function, 
!    SIAM Journal on Numerical Analysis, B, 
!    Volume 1, pages 86-96, 1964.
!
!  Modified:
!
!    30 March 1999
!
!  Author:
!
!    Alan Miller
!    CSIRO Division of Mathematics and Statistics
!
!  Parameters:
!
!    Input, real Z, the argument of the Gamma function.
!
!    Output, integer IER, error flag.
!    0, no error occurred.
!    1, Z is less than or equal to 0.
!
!    Output, real LNGAMMA, the logarithm of the gamma function of Z.
!
  double precision, parameter :: lnsqrt2pi = 0.9189385332046727
!
  double precision, parameter, dimension ( 9 ) :: a = (/ &
            0.9999999999995183, &
          676.5203681218835, &
       - 1259.139216722289, &
          771.3234287757674, &
        - 176.6150291498386, &
           12.50734324009056, &
          - 0.1385710331296526, &
            0.9934937113930748E-05, &
            0.1659470187408462E-06 /)
  integer ier
  integer j
  double precision lngamma
  double precision tmp
  double precision z
!
  if ( z <= 0.0 ) then
    ier = 1
    lngamma = 0.0
    write ( *, * ) ' '
    write ( *, * ) 'LNGAMMA - Fatal error!'
    write ( *, * ) '  Z <= 0.'
    stop
  end if

  ier = 0

  lngamma = 0.0
  tmp = z + 7.0
  do j = 9, 2, -1
    lngamma = lngamma + a(j) / tmp
    tmp = tmp - 1.0
  end do

  lngamma = lngamma + a(1)
  lngamma = log ( lngamma ) + lnsqrt2pi - ( z + 6.5 ) + ( z - 0.5 ) * &
    log ( z + 6.5 )

  return
end
subroutine normal_01_sample ( iseed, x )
!
!*******************************************************************************
!
!! NORMAL_01_SAMPLE samples the standard Normal PDF.
!
!
!  Discussion:
!
!    The standard normal distribution has mean 0 and standard
!    deviation 1.
!
!  Method:
!
!    The Box-Muller method is used.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  double precision, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  integer iseed
  integer iset
  double precision uniform_01_sample
  double precision v1
  double precision v2
  double precision x
  double precision xsave
!
  save iset
  save xsave

  data iset / 0 /
  data xsave / 0.0 /
!
  if ( iset == 0 ) then

    v1 = uniform_01_sample ( iseed )

    if ( v1 <= 0.0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'NORMAL_01_SAMPLE - Fatal error!'
      write ( *, * ) '  V1 <= 0.'
      write ( *, * ) '  V1 = ', v1
      write ( *, * ) '  ISEED = ', iseed
      stop
    end if

    v2 = uniform_01_sample ( iseed )

    if ( v2 <= 0.0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'NORMAL_01_SAMPLE - Fatal error!'
      write ( *, * ) '  V2 <= 0.'
      write ( *, * ) '  V2 = ', v2
      write ( *, * ) '  ISEED = ', iseed
      stop
    end if

    x = sqrt ( - 2.0 * log ( v1 ) ) * cos ( 2.0 * PI * v2 )

    xsave = sqrt ( - 2.0 * log ( v1 ) ) * sin ( 2.0 * PI * v2 )

    iset = 1

  else

    x = xsave
    iset = 0

  end if

  return
end
subroutine normp ( z, p, q, pdf )
!
!*******************************************************************************
!
!! NORMP computes the cumulative density of the standard normal distribution.
!
!
!  Reference:
!
!    Algorithm 5666 for the error function,
!    Hart, Cheney, Lawson, Machly, Mesztenyi, Rice, Thacher and Witzgall, 
!    Computer Approximations, 
!    Wiley, 1968.
!
!  Author:
!
!    Alan Miller
!
!  Modified:
!
!    27 March 1999
!
!  Parameters:
!
!    Input, real Z, divides the real line into two semi-infinite
!    intervals, over each of which the standard normal distribution
!    is to be integrated.
!
!    Output, real P, Q, the integrals of the standard normal
!    distribution over the intervals ( - Infinity, Z] and 
!    [Z, + Infinity ), respectively.
!
!    Output, real PDF, the value of the standard normal distribution at Z.
!
  double precision, parameter :: cutoff = 7.071
  double precision, parameter :: p0 = 220.2068679123761
  double precision, parameter :: p1 = 221.2135961699311
  double precision, parameter :: p2 = 112.0792914978709
  double precision, parameter :: p3 = 33.91286607838300
  double precision, parameter :: p4 = 6.373962203531650
  double precision, parameter :: p5 = 0.7003830644436881
  double precision, parameter :: p6 = 0.03526249659989109
  double precision, parameter :: q0 = 440.4137358247522
  double precision, parameter :: q1 = 793.8265125199484
  double precision, parameter :: q2 = 637.3336333788311
  double precision, parameter :: q3 = 296.5642487796737
  double precision, parameter :: q4 = 86.78073220294608
  double precision, parameter :: q5 = 16.06417757920695
  double precision, parameter :: q6 = 1.755667163182642
  double precision, parameter :: q7 = 0.08838834764831844
  double precision, parameter :: root2pi = 2.506628274631001
!
  double precision expntl
  double precision p
  double precision pdf
  double precision q
  double precision z
  double precision zabs
!
  zabs = abs ( z )
!
!  |Z| > 37.
!
  if ( zabs > 37.0 ) then

    pdf = 0.0
    p = 0.0
!
!  |Z| <= 37.
!
  else

    expntl = exp ( - 0.5 * zabs**2 )
    pdf = expntl / root2pi
!
!  |Z| < CUTOFF = 10 / sqrt(2).
!
    if ( zabs < cutoff ) then

      p = expntl * (((((( &
             p6   * zabs &
           + p5 ) * zabs &
           + p4 ) * zabs &
           + p3 ) * zabs &
           + p2 ) * zabs &
           + p1 ) * zabs &
           + p0 ) / ((((((( &
             q7   * zabs &
           + q6 ) * zabs &
           + q5 ) * zabs &
           + q4 ) * zabs &
           + q3 ) * zabs &
           + q2 ) * zabs &
           + q1 ) * zabs &
           + q0 )
!
!  |Z| >= CUTOFF.
!
    else

      p = pdf / ( &
           zabs + 1.0 / ( &
           zabs + 2.0 / ( &
           zabs + 3.0 / ( &
           zabs + 4.0 / ( &
           zabs + 0.65 )))))

    end if

  end if

  if ( z < 0.0 ) then
    q = 1.0 - p
  else
    q = p
    p = 1.0 - q
  end if

  return
end
subroutine nprob ( z, p, q, pdf )
!
!*******************************************************************************
!
!! NPROB computes the cumulative density of the standard normal distribution.
!
!
!  Reference:
!
!    A G Adams,
!    Algorithm 39, Areas Under the Normal Curve,
!    Computer Journal, 
!    Volume 12, pages 197-198, 1969.
!
!  Modified:
!
!    31 March 1999
!
!  Parameters:
!
!    Input, real Z, divides the real line into two semi-infinite
!    intervals, over each of which the standard normal distribution
!    is to be integrated.
!
!    Output, real P, Q, the integrals of the standard normal
!    distribution over the intervals ( - Infinity, Z] and 
!    [Z, + Infinity ), respectively.
!
!    Output, real PDF, the value of the standard normal distribution at Z.
!
  double precision, parameter :: a0 = 0.5
  double precision, parameter :: a1 = 0.398942280444
  double precision, parameter :: a2 = 0.399903438504
  double precision, parameter :: a3 = 5.75885480458
  double precision, parameter :: a4 = 29.8213557808
  double precision, parameter :: a5 = 2.62433121679
  double precision, parameter :: a6 = 48.6959930692
  double precision, parameter :: a7 = 5.92885724438
  double precision, parameter :: b0 = 0.398942280385
  double precision, parameter :: b1 = 0.000000038052
  double precision, parameter :: b2 = 1.00000615302
  double precision, parameter :: b3 = 0.000398064794
  double precision, parameter :: b4 = 1.98615381364
  double precision, parameter :: b5 = 0.151679116635 
  double precision, parameter :: b6 = 5.29330324926
  double precision, parameter :: b7 = 4.8385912808
  double precision, parameter :: b8 = 15.1508972451
  double precision, parameter :: b9 = 0.742380924027
  double precision, parameter :: b10 = 30.789933034
  double precision, parameter :: b11 = 3.99019417011
!
  double precision p
  double precision pdf
  double precision q
  double precision y
  double precision z
  double precision zabs
!
  zabs = abs ( z )
!
!  |Z| between 0 and 1.28
!
  if ( abs ( z ) <= 1.28 ) then

    y = a0 * z**2
    pdf = exp ( - y ) * b0

    q = a0 - zabs * ( a1 - a2 * y &
         / ( y + a3 - a4 &
         / ( y + a5 + a6 &
         / ( y + a7 ))))
!
!  |Z| between 1.28 and 12.7
!
  else if ( abs ( z ) <= 12.7 ) then

    y = a0 * z**2
    pdf = exp ( - y ) * b0

    q = pdf &
         / ( zabs - b1 + b2 &
         / ( zabs + b3 + b4 &
         / ( zabs - b5 + b6 &
         / ( zabs + b7 - b8 &
         / ( zabs + b9 + b10 &
         / ( zabs + b11 ))))))
!
!  Z far out in tail.
!
  else

    q = 0.0
    pdf = 0.0

  end if

  if ( z < 0.0 ) then
    p = q
    q = 1.0 - p
  else
    p = 1.0 - q
  end if

  return
end
function ppchi2 ( p, v, ifault )
!
!*******************************************************************************
!
!! PPCHI2 evaluates the percentage points of the chi-squared PDF.
!
!
!  Auxiliary routines:
!
!     PPND = AS 111 or AS 241;
!     GAMMAD = AS 239.
!
!  Reference:
!
!    Best and Roberts,
!    The Percentage Points of the Chi-Squared Distribution,
!    Algorithm AS 91,  
!    Applied Statistics,
!    Volume 24, Number ?, pages 385-390, 1975.
!
!  Modified:
!
!    30 March 1999
!
!  Parameters:
!
!    Input, real P, a value of the chi-squared cumulative probability
!    density function.
!    0.000002 <= P <= 0.999998.
!
!    Input, real V, the parameter of the chi-squared probability density
!    function.  V > 0.
!
!    Output, integer IFAULT, error flag.
!    0, no error detected.
!    1, P < PMIN or P > PMAX.
!    2, V <= 0.0.
!    3, an error occurred in the incomplete Gamma function routine.
!    4, the maximum number of iterations were taken without convergence.
!    5, an error occurred in the log Gamma routine.
!
!    Output, real PPCHI2, the value of the chi-squared random deviate
!    with the property that the probability that a chi-squared random
!    deviate with parameter V is less than or equal to PPCHI2 is P.
!

  double precision, parameter :: aa = 0.6931471806
  double precision, parameter :: c1 = 0.01
  double precision, parameter :: c2 = 0.222222
  double precision, parameter :: c3 = 0.32
  double precision, parameter :: c4 = 0.4
  double precision, parameter :: c5 = 1.24
  double precision, parameter :: c6 = 2.2
  double precision, parameter :: c7 = 4.67
  double precision, parameter :: c8 = 6.66
  double precision, parameter :: c9 = 6.73
  double precision, parameter :: c10 = 13.32
  double precision, parameter :: c11 = 60.0
  double precision, parameter :: c12 = 70.0
  double precision, parameter :: c13 = 84.0
  double precision, parameter :: c14 = 105.0
  double precision, parameter :: c15 = 120.0
  double precision, parameter :: c16 = 127.0
  double precision, parameter :: c17 = 140.0
  double precision, parameter :: c18 = 175.0
  double precision, parameter :: c19 = 210.0
  double precision, parameter :: c20 = 252.0
  double precision, parameter :: c21 = 264.0
  double precision, parameter :: c22 = 294.0
  double precision, parameter :: c23 = 346.0
  double precision, parameter :: c24 = 420.0
  double precision, parameter :: c25 = 462.0
  double precision, parameter :: c26 = 606.0
  double precision, parameter :: c27 = 672.0
  double precision, parameter :: c28 = 707.0
  double precision, parameter :: c29 = 735.0
  double precision, parameter :: c30 = 889.0
  double precision, parameter :: c31 = 932.0
  double precision, parameter :: c32 = 966.0
  double precision, parameter :: c33 = 1141.0
  double precision, parameter :: c34 = 1182.0
  double precision, parameter :: c35 = 1278.0
  double precision, parameter :: c36 = 1740.0
  double precision, parameter :: c37 = 2520.0
  double precision, parameter :: c38 = 5040.0
  double precision, parameter :: e = 0.0000005
  integer, parameter :: maxit = 20
  double precision, parameter :: pmax = 0.999998
  double precision, parameter :: pmin = 0.000002
!
  double precision a
  double precision alngam
  double precision b
  double precision c
  double precision ch
  double precision g
  double precision gammad
  integer i
  integer ifault
  integer ifault2
  double precision p
  double precision p1
  double precision p2
  double precision ppchi2
  double precision ppnd
  double precision q
  double precision s1
  double precision s2
  double precision s3
  double precision s4
  double precision s5
  double precision s6
  double precision t
  double precision v
  double precision x
  double precision xx
!
  ifault2 = 0
!
!  Check the input.
!
  if ( p < pmin ) then
    ifault = 1
    ppchi2 = - 1.0
    write ( *, * ) ' '
    write ( *, * ) 'PPCHI2 - Fatal error!'
    write ( *, * ) '  P < PMIN.'
    stop
  end if

  if ( p > pmax ) then
    ifault = 1
    ppchi2 = - 1.0
    write ( *, * ) ' '
    write ( *, * ) 'PPCHI2 - Fatal error!'
    write ( *, * ) '  P > PMAX.'
    stop
  end if

  if ( v <= 0.0 ) then
    ifault = 2
    ppchi2 = - 1.0
    write ( *, * ) ' '
    write ( *, * ) 'PPCHI2 - Fatal error!'
    write ( *, * ) '  V <= 0.0.'
    stop
  end if

  ifault = 0
  xx = 0.5 * v
  c = xx - 1.0
!
!  Compute Log ( Gamma ( V/2 ) ).
!
  g = alngam ( v / 2.0, ifault )

  if ( ifault /= 0 ) then
    ifault = 5
    write ( *, * ) ' '
    write ( *, * ) 'PPCHI2 - Fatal error!'
    write ( *, * ) '  ALNGAM returns error code.'
    stop
  end if
!
!  Starting approximation for small chi-squared.
!
  if ( v < - c5 * log ( p ) ) then

    ch = ( p * xx * exp ( g + xx * aa ) )**( 1.0 / xx )

    if ( ch < e ) then
      ifault = 0
      ppchi2 = ch
      return
    end if
!
!  Starting approximation for V less than or equal to 0.32.
!
  else if ( v <= c3 ) then

    ch = c4
    a = log ( 1.0 - p )

10      continue

    q = ch
    p1 = 1.0 + ch * ( c7 + ch )
    p2 = ch * ( c9 + ch * ( c8 + ch ) )

    t = - 0.5 + ( c7 + 2.0 * ch ) / p1 - ( c9 + ch * ( c10 + 3.0 * ch ) ) / p2

    ch = ch - ( 1.0 - exp ( a + g + 0.5 * ch + c * aa ) * p2 / p1 ) / t

    if ( abs ( q / ch - 1.0 ) > c1 ) then
      go to 10
    end if
!
!  Call to algorithm AS 111.
!  Note that P has been tested above.
!  AS 241 could be used as an alternative.
!
  else

    x = ppnd ( p, ifault2 )
!
!  Starting approximation using Wilson and Hilferty estimate.
!
    p1 = c2 / v
    ch = v * ( x * sqrt ( p1 ) + 1.0 - p1 )**3
!
!  Starting approximation for P tending to 1.
!
    if ( ch > c6 * v + 6.0 ) then
      ch = - 2.0 * ( log ( 1.0 - p ) - c * log ( 0.5 * ch ) + g )
    end if

  end if
!
!  Call to algorithm AS 239 and calculation of seven term Taylor series.
!
  do i = 1, maxit

    q = ch
    p1 = 0.5 * ch
    p2 = p - gammad ( p1, xx, ifault2 )

    if ( ifault2 /= 0 ) then
      ppchi2 = - 1.0
      ifault = 3
      write ( *, * ) ' '
      write ( *, * ) 'PPCHI2 - Fatal error!'
      write ( *, * ) '  GAMMAD returns error code.'
      stop
    end if

    t = p2 * exp ( xx * aa + g + p1 - c * log ( ch ) )
    b = t / ch
    a = 0.5 * t - b * c

    s1 = &
           ( c19 + a &
         * ( c17 + a &
         * ( c14 + a &
         * ( c13 + a &
         * ( c12 + a &
         *   c11 ) ) ) ) ) / c24

    s2 = &
           ( c24 + a &
         * ( c29 + a &
         * ( c32 + a &
         * ( c33 + a &
         *   c35 ) ) ) ) / c37

    s3 = &
           ( c19 + a &
         * ( c25 + a &
         * ( c28 + a * &
             c31 ) ) ) / c37

    s4 = &
           ( c20 + a &
         * ( c27 + a &
         *   c34 ) + c &
         * ( c22 + a &
         * ( c30 + a &
         *   c36 ) ) ) / c38

    s5 = ( c13 + c21 * a + c * ( c18 + c26 * a ) ) / c37

    s6 = ( c15 + c * ( c23 + c16 * c ) ) / c38

    ch = ch + t * ( 1.0 + 0.5 * t * s1 - b * c &
         * ( s1 - b &
         * ( s2 - b &
         * ( s3 - b &
         * ( s4 - b &
         * ( s5 - b &
         *   s6 ) ) ) ) ) )

    if ( abs ( q / ch - 1.0 ) > e ) then
      ifault = 0
      ppchi2 = ch
      return
    end if

  end do

  ifault = 4
  ppchi2 = ch
  write ( *, * ) ' '
  write ( *, * ) 'PPCHI2 - Warning!'
  write ( *, * ) '  Convergence not reached.'

  return
end
function ppnd ( p, ifault )
!
!*******************************************************************************
!
!! PPND produces the normal deviate value corresponding to lower tail area = P.
!
!
!  Reference:
!
!    J Beasley and S Springer,
!    The Percentage Points of the Normal Distribution,
!    Algorithm AS 111, 
!    Applied Statistics,
!    Volume 26, Number ?, pages 118-121, 1977.
!
!  Modified:
!
!    28 March 1999
!
!  Parameters:
!
!    Input, real P, the value of the cumulative probability densitity function.
!    0 < P < 1.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    1, P <= 0 or P >= 1.  PPND is returned as 0.
!
!    Output, real PPND, the normal deviate value with the property that
!    the probability of a standard normal deviate being less than or
!    equal to PPND is P.
!
  double precision, parameter :: a0 = 2.50662823884
  double precision, parameter :: a1 = -18.61500062529
  double precision, parameter :: a2 = 41.39119773534
  double precision, parameter :: a3 = -25.44106049637
  double precision, parameter :: b1 = -8.47351093090
  double precision, parameter :: b2 = 23.08336743743
  double precision, parameter :: b3 = -21.06224101826
  double precision, parameter :: b4 = 3.13082909833
  double precision, parameter :: c0 = -2.78718931138
  double precision, parameter :: c1 = -2.29796479134
  double precision, parameter :: c2 = 4.85014127135
  double precision, parameter :: c3 = 2.32121276858
  double precision, parameter :: d1 = 3.54388924762
  double precision, parameter :: d2 = 1.63706781897
  double precision, parameter :: split = 0.42
!
  integer ifault
  double precision p
  double precision ppnd
  double precision r
!
  ifault = 0
!
!  0.08 < P < 0.92
!
  if ( abs ( p - 0.5 ) <= split ) then

    r = ( p - 0.5 )**2

    ppnd = ( p - 0.5 ) * ( ( ( &
           a3   * r &
         + a2 ) * r &
         + a1 ) * r &
         + a0 ) / ( ( ( ( &
           b4   * r &
         + b3 ) * r &
         + b2 ) * r &
         + b1 ) * r &
         + 1.0 )
!
!  P < 0.08 or P > 0.92, 
!  R = min ( P, 1-P )
!
  else if ( p > 0.0 .and. p < 1.0 ) then

    if ( p > 0.5 ) then
      r = sqrt ( - log ( 1.0 - p ) )
    else
      r = sqrt ( - log ( p ) )
    end if

    ppnd = ( ( ( &
           c3   * r &
         + c2 ) * r &
         + c1 ) * r &
         + c0 ) / ( ( &
           d2   * r &
         + d1 ) * r &
         + 1.0 )

    if ( p < 0.5 ) then
      ppnd = - ppnd
    end if
!
!  P <= 0.0 or P >= 1.0
!
  else

    ifault = 1
    ppnd = 0.0

    write ( *, * ) ' '
    write ( *, * ) 'PPND - Warning!'
    write ( *, * ) '  P <= 0 or P >=1.'
    write ( *, * ) '  PPND value would be infinite.'
    
  end if

  return
end
function ppnd16 ( p, ifault )
!
!*******************************************************************************
!
!! PPND16 produces the normal deviate value corresponding to lower tail area = P.
!
!
!  Accuracy:
!
!    The result is accurate to about 1 part in 10**16.
!
!  Reference:
!
!    Michael Wichura,
!    The Percentage Points of the Normal Distribution,
!    Algorithm AS 241,
!    Applied Statistics,
!    Volume 37, Number 3, pages 477-484, 1988.
!
!  Modified:
!
!    28 March 1999
!
!  Parameters:
!
!    Input, double precision P, the value of the cumulative probability 
!    densitity function.  0 < P < 1.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    1, P <= 0 or P >= 1.
!
!    Output, double precision PPND16, the normal deviate value with the 
!    property that the probability of a standard normal deviate being 
!    less than or equal to PPND16 is P.
!
  double precision, parameter :: a0 = 3.3871328727963666080d0
  double precision, parameter :: a1 = 1.3314166789178437745d+2
  double precision, parameter :: a2 = 1.9715909503065514427d+3
  double precision, parameter :: a3 = 1.3731693765509461125d+4
  double precision, parameter :: a4 = 4.5921953931549871457d+4
  double precision, parameter :: a5 = 6.7265770927008700853d+4
  double precision, parameter :: a6 = 3.3430575583588128105d+4
  double precision, parameter :: a7 = 2.5090809287301226727d+3
  double precision, parameter :: b1 = 4.2313330701600911252d+1
  double precision, parameter :: b2 = 6.8718700749205790830d+2
  double precision, parameter :: b3 = 5.3941960214247511077d+3
  double precision, parameter :: b4 = 2.1213794301586595867d+4
  double precision, parameter :: b5 = 3.9307895800092710610d+4
  double precision, parameter :: b6 = 2.8729085735721942674d+4
  double precision, parameter :: b7 = 5.2264952788528545610d+3
  double precision, parameter :: c0 = 1.42343711074968357734d0
  double precision, parameter :: c1 = 4.63033784615654529590d0
  double precision, parameter :: c2 = 5.76949722146069140550d0
  double precision, parameter :: c3 = 3.64784832476320460504d0
  double precision, parameter :: c4 = 1.27045825245236838258d0
  double precision, parameter :: c5 = 2.41780725177450611770d-1
  double precision, parameter :: c6 = 2.27238449892691845833d-2
  double precision, parameter :: c7 = 7.74545014278341407640d-4
  double precision, parameter :: const1 = 0.180625d0
  double precision, parameter :: const2 = 1.6d0
  double precision, parameter :: d1 = 2.05319162663775882187d0
  double precision, parameter :: d2 = 1.67638483018380384940d0
  double precision, parameter :: d3 = 6.89767334985100004550d-1
  double precision, parameter :: d4 = 1.48103976427480074590d-1
  double precision, parameter :: d5 = 1.51986665636164571966d-2
  double precision, parameter :: d6 = 5.47593808499534494600d-4
  double precision, parameter :: d7 = 1.05075007164441684324d-9
  double precision, parameter :: e0 = 6.65790464350110377720d0
  double precision, parameter :: e1 = 5.46378491116411436990d0
  double precision, parameter :: e2 = 1.78482653991729133580d0
  double precision, parameter :: e3 = 2.96560571828504891230d-1
  double precision, parameter :: e4 = 2.65321895265761230930d-2
  double precision, parameter :: e5 = 1.24266094738807843860d-3
  double precision, parameter :: e6 = 2.71155556874348757815d-5
  double precision, parameter :: e7 = 2.01033439929228813265d-7
  double precision, parameter :: f1 = 5.99832206555887937690d-1
  double precision, parameter :: f2 = 1.36929880922735805310d-1
  double precision, parameter :: f3 = 1.48753612908506148525d-2
  double precision, parameter :: f4 = 7.86869131145613259100d-4
  double precision, parameter :: f5 = 1.84631831751005468180d-5
  double precision, parameter :: f6 = 1.42151175831644588870d-7
  double precision, parameter :: f7 = 2.04426310338993978564d-15
  double precision, parameter :: split1 = 0.425d0
  double precision, parameter :: split2 = 5.d0
!
  integer ifault
  double precision p
  double precision ppnd16
  double precision q
  double precision r
!
  ifault = 0
  q = p - 0.5

  if ( abs ( q ) <= split1 ) then

    r = const1 - q**2

    ppnd16 = q * ((((((( &
           a7   * r &
         + a6 ) * r &
         + a5 ) * r &
         + a4 ) * r &
         + a3 ) * r &
         + a2 ) * r &
         + a1 ) * r &
         + a0 ) / ((((((( &
           b7   * r &
         + b6 ) * r &
         + b5 ) * r &
         + b4 ) * r &
         + b3 ) * r &
         + b2 ) * r &
         + b1 ) * r &
         + 1.0 )

  else

    if ( q < 0.0 ) then
      r = p
    else
      r = 1.0 - p
    end if

    if ( r <= 0.0 ) then
      ifault = 1
      ppnd16 = 0.0
      write ( *, * ) ' '
      write ( *, * ) 'PPND16 - Warning!'
      write ( *, * ) '  P <= 0 or P >= 1.'
      write ( *, * ) '  PPND16 value would be infinite.'
      return
    end if

    r = sqrt ( - log ( r ) )

    if ( r <= split2 ) then

      r = r - const2

      ppnd16 = ((((((( &
             c7   * r &
           + c6 ) * r &
           + c5 ) * r &
           + c4 ) * r &
           + c3 ) * r &
           + c2 ) * r &
           + c1 ) * r &
           + c0 ) / ((((((( &
             d7   * r &
           + d6 ) * r &
           + d5 ) * r &
           + d4 ) * r &
           + d3 ) * r &
           + d2 ) * r &
           + d1 ) * r &
           + 1.0 )

    else

      r = r - split2

      ppnd16 = ((((((( &
             e7   * r &
           + e6 ) * r &
           + e5 ) * r &
           + e4 ) * r &
           + e3 ) * r &
           + e2 ) * r &
           + e1 ) * r &
           + e0 ) / ((((((( &
             f7   * r &
           + f6 ) * r &
           + f5 ) * r &
           + f4 ) * r &
           + f3 ) * r &
           + f2 ) * r &
           + f1 ) * r &
           + 1.0 )

    end if

    if ( q < 0.0 ) then
      ppnd16 = - ppnd16
    end if

  end if

  return
end
function ppnd7 ( p, ifault )
!
!*******************************************************************************
!
!! PPND7 produces the normal deviate value corresponding to lower tail area = P.
!
!
!  Accuracy:
!
!    The result is accurate to about 1 part in 10**7.
!
!  Reference:
!
!    Michael Wichura,
!    The Percentage Points of the Normal Distribution,
!    Algorithm AS 241,
!    Applied Statistics,
!    Volume 37, Number 3, pages 477-484, 1988.
!
!  Modified:
!
!    31 March 1999
!
!  Parameters:
!
!    Input, real P, the value of the cumulative probability densitity function.
!    0 < P < 1.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    1, P <= 0 or P >= 1.
!
!    Output, real PPND7, the normal deviate value with the property that
!    the probability of a standard normal deviate being less than or
!    equal to PPND7 is P.
!
  double precision, parameter :: a0 = 3.3871327179
  double precision, parameter :: a1 = 50.434271938
  double precision, parameter :: a2 = 159.29113202
  double precision, parameter :: a3 = 59.109374720
  double precision, parameter :: b1 = 17.895169469
  double precision, parameter :: b2 = 78.757757664
  double precision, parameter :: b3 = 67.187563600
  double precision, parameter :: c0 = 1.4234372777
  double precision, parameter :: c1 = 2.7568153900
  double precision, parameter :: c2 = 1.3067284816
  double precision, parameter :: c3 = 0.17023821103
  double precision, parameter :: const1 = 0.180625
  double precision, parameter :: const2 = 1.6
  double precision, parameter :: d1 = 0.73700164250
  double precision, parameter :: d2 = 0.12021132975
  double precision, parameter :: e0 = 6.6579051150
  double precision, parameter :: e1 = 3.0812263860
  double precision, parameter :: e2 = 0.42868294337
  double precision, parameter :: e3 = 0.017337203997
  double precision, parameter :: f1 = 0.24197894225
  double precision, parameter :: f2 = 0.012258202635
  double precision, parameter :: split1 = 0.425
  double precision, parameter :: split2 = 5.0
!
  integer ifault
  double precision p
  double precision ppnd7
  double precision q
  double precision r
!
  ifault = 0
  q = p - 0.5

  if ( abs ( q ) <= split1 ) then

    r = const1 - q**2

    ppnd7 = q * ((( &
           a3   * r &
         + a2 ) * r &
         + a1 ) * r &
         + a0 ) / ((( &
           b3   * r &
         + b2 ) * r &
         + b1 ) * r &
         + 1.0 )

  else

    if ( q < 0.0 ) then
      r = p
    else
      r = 1.0 - p
    end if

    if ( r <= 0.0 ) then
      ifault = 1
      ppnd7 = 0.0
      write ( *, * ) ' '
      write ( *, * ) 'PPND7 - Fatal error!'
      write ( *, * ) '  P <= 0 or P >= 1.'
      write ( *, * ) '  PPND7 value would be infinite.'
      return
    end if

    r = sqrt ( - log ( r ) )

    if ( r <= split2 ) then

      r = r - const2

      ppnd7 = ((( &
           c3   * r &
         + c2 ) * r &
         + c1 ) * r &
         + c0 ) / (( &
           d2   * r &
         + d1 ) * r &
         + 1.0 )

    else

      r = r - split2

      ppnd7 = ((( &
           e3   * r &
         + e2 ) * r &
         + e1 ) * r &
         + e0 ) / (( &
           f2   * r &
         + f1 ) * r &
         + 1.0 )

    end if

    if ( q < 0.0 ) then
      ppnd7 = - ppnd7
    end if

  end if

  return
end
function psi ( xx )
!
!*******************************************************************************
!
!! PSI evaluates the psi or digamma function, d/dx ln(gamma(x)).
!
!
!  Discussion:
!
!    The main computation involves evaluation of rational Chebyshev
!    approximations published in Math. Comp. 27, 123-127(1973) by
!    Cody, Strecok and Thacher.
!
!    PSI was written at Argonne National Laboratory for the FUNPACK
!    package of special function subroutines.  PSI was modified by
!    A. H. Morris (NSWC).
!
!  Parameters:
!
!    Input, real XX, the argument of the psi function.
!
!    Output, real PSI, the value of the psi function.  PSI is assigned
!    the value 0 when the psi function is undefined.
!
  double precision, parameter :: dx0 = 1.461632144968362341262659542325721325
  double precision, parameter :: piov4 = 0.785398163397448
!
  double precision aug
  double precision den
  double precision eps
  integer i
  integer m
  integer n
  integer nq
  double precision, parameter, dimension ( 7 ) :: p1 = (/&
    0.895385022981970e-02, &
    0.477762828042627e+01, &
    0.142441585084029e+03, &
    0.118645200713425e+04, &
    0.363351846806499e+04, &
    0.413810161269013e+04, &
    0.130560269827897e+04 /)
  double precision, parameter, dimension ( 4 ) :: p2 = (/ &
    -0.212940445131011e+01, &
    -0.701677227766759e+01, &
    -0.448616543918019e+01, &
    -0.648157123766197e+00 /)
  double precision psi
  double precision, parameter, dimension ( 6 ) :: q1 = (/&
    0.448452573429826e+02, &
    0.520752771467162e+03, &
    0.221000799247830e+04, &
    0.364127349079381e+04, &
    0.190831076596300e+04, &
    0.691091682714533e-05 /)
  double precision, parameter, dimension ( 4 ) :: q2 = (/ &
    0.322703493791143e+02, &
    0.892920700481861e+02, &
    0.546117738103215e+02, &
    0.777788548522962e+01 /)
  double precision sgn
  double precision temp
  double precision upper
  double precision w
  double precision x
  double precision xmax1
  double precision xmx0
  double precision xsmall
  double precision xx
  double precision z
!
!  XMAX1 is the largest positive floating point constant with entirely
!  integer representation.  It is also used as negative of lower bound
!  on acceptable negative arguments and as the positive argument beyond which
!  psi may be represented as LOG(X).
!
  xmax1 = dble ( 2147483647 )

  eps = 1.0
10    continue
  eps = eps / 2.0
  temp = eps + 1.0
  if ( temp > 1.0 ) then
    go to 10
  end if
  eps = 2.0 * eps

  xmax1 = min ( xmax1, 1.0 / eps )
!
!  XSMALL is the absolute argument below which PI*COTAN(PI*X)
!  may be represented by 1/X.
!
  xsmall = 1.0e-9

  x = xx
  aug = 0.0

  if ( x == 0.0 ) then
    psi = 0.0
    return
  end if
!
!  X < 0.5,  Use reflection formula PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
!
  if ( x < 0.5 ) then
!
!  0 < ABS(X) <= XSMALL.  Use 1/X as a substitute for PI*COTAN(PI*X)
!
    if ( abs ( x ) <= xsmall ) then
      aug = - 1.0 / x
      go to 40
    end if
!
!  Reduction of argument for cotan
!
    w = - x
    sgn = piov4

    if ( w <= 0.0 ) then
      w = - w
      sgn = -sgn
    end if
!
!  Make an error exit if X .LE. -XMAX1
!
    if ( w >= xmax1 ) then
      psi = 0.0
      return
    end if

    nq = int ( w )
    w = w - dble ( nq )
    nq = int ( w * 4.0 )
    w = 4.0 * ( w - dble ( nq ) * 0.25 )
!
!  W is now related to the fractional part of  4.0 * X.
!  Adjust argument to correspond to values in first
!  quadrant and determine sign
!
    n = nq / 2
    if ( ( n + n ) /= nq ) then
      w = 1.0 - w
    end if

    z = piov4 * w
    m = n / 2

    if ( ( m + m ) /= n ) then
      sgn = - sgn
    end if
!
!  Determine final value for  -PI*COTAN(PI*X)
!
    n = ( nq + 1 ) / 2
    m = n / 2
    m = m + m

    if ( m == n ) then

      if ( z == 0.0 ) then
        psi = 0.0
       return
    end  if

      aug = 4.0 * sgn * ( cos(z) / sin(z) )

    else

      aug = 4.0 * sgn * ( sin(z) / cos(z) )

    end if

   40   continue

    x = 1.0 - x

  end if
!
!  0.5 <= X <= 3.0
!
  if ( x <= 3.0 ) then

    den = x
    upper = p1(1) * x

    do i = 1, 5
      den = ( den + q1(i) ) * x
      upper = ( upper + p1(i+1) ) * x
    end do

    den = ( upper + p1(7) ) / ( den + q1(6) )
    xmx0 = dble ( x ) - dx0
    psi = den * xmx0 + aug
!
!  3.0 < X < XMAX1
!
  else if ( x < xmax1 ) then

    w = 1.0 / x**2
    den = w
    upper = p2(1) * w

    do i = 1, 3
      den = ( den + q2(i) ) * w
      upper = ( upper + p2(i+1) ) * w
    end do

    aug = upper / ( den + q2(4) ) - 0.5 / x + aug
    psi = aug + log ( x )
!
!  XMAX1 <= X
!
  else

    psi = aug + log ( x )

  end if

  return
end
subroutine r_random ( r, rlo, rhi, iseed )
!
!*******************************************************************************
!
!! R_RANDOM returns a random real in a given range.
!
!
!  Modified:
!
!    26 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real R, the randomly chosen value.
!
!    Input, real RLO, RHI, the minimum and maximum values.
!
!    Input/output, integer ISEED, a seed for the random number generator.
!
  integer iseed
  double precision r
  double precision rhi
  double precision rlo
  double precision t
  double precision uniform_01_sample
!
!  Pick a random number in (0,1).
!
  t = uniform_01_sample ( iseed )
!
!  Set R.
!
  r = ( 1.0 - t ) * rlo + t * rhi

  return
end
subroutine rcol_mean ( lda, m, n, a, mean )
!
!*******************************************************************************
!
!! RCOL_MEAN returns the means of columns of a real array.
!
!
!  Example:
!
!    A =
!      1  2  3
!      2  6  7
!
!    MEAN =
!      1.5  4.0  5.0
!
!  Modified:
!
!    29 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array, which should
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real A(LDA,N), the array to be examined.
!
!    Output, real MEAN(N), the means, or averages, of the columns.
!
  integer lda
  integer m
  integer n
!
  double precision a(lda,n)
  integer i
  integer j
  double precision mean(n)
!
  do j = 1, n

    mean(j) = 0.0
    do i = 1, m
      mean(j) = mean(j) + a(i,j)
    end do

    mean(j) = mean(j) / dble ( m )

  end do

  return
end
subroutine rcol_variance ( lda, m, n, a, variance )
!
!*******************************************************************************
!
!! RCOL_VARIANCE returns the variances of the columns of a real array.
!
!
!  Modified:
!
!    29 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A, which should
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns in the array.
!
!    Input, real A(LDA,N), the array whose variances are desired.
!
!    Output, real VARIANCE(N), the variances of the rows.
!
  integer lda
  integer m
  integer n
!
  double precision a(lda,n)
  integer i
  integer j
  double precision mean
  double precision variance(n)
!
  do j = 1, n

    mean = 0.0
    do i = 1, m
      mean = mean + a(i,j)
    end do
    mean = mean / dble ( m )

    variance(j) = 0.0
    do i = 1, m
      variance(j) = variance(j) + ( a(i,j) - mean )**2
    end do

    if ( m > 1 ) then
      variance(j) = variance(j) / dble ( m - 1 )
    else
      variance(j) = 0.0
    end if

  end do

  return
end
subroutine rvec_sum ( n, a, sum )
!
!*******************************************************************************
!
!! RVEC_SUM sums the entries of a real vector.
!
!
!  Modified:
!
!    04 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real A(N), the vector whose sum is desired.
!
!    Output, real SUM, the sum of the vector entries.
!
  integer n
!
  double precision a(n)
  integer i
  double precision sum
!
  sum = 0.0
  do i = 1, n
    sum = sum + a(i)
  end do

  return
end
function trigamma ( x )
!
!*******************************************************************************
!
!! TRIGAMMA calculates trigamma(x) = d**2 log(Gamma(x)) / dx**2.
!
!
!  Reference:
!
!    B Schneider,
!    Trigamma Function,
!    Algorithm AS 121,
!    Applied Statistics, 
!    Volume 27, Number 1, page 97-99, 1978.
!
!  Modified:
!
!    03 January 2000
!
!  Parameters:
!
!    Input, real X, the argument of the trigamma function.
!    0 < X.
!
!    Output, real TRIGAMMA, the value of the trigamma function at X.
!
  double precision, parameter :: a = 0.0001
  double precision, parameter :: b = 5.0
  double precision, parameter :: b2 =   1.0 / 6.0
  double precision, parameter :: b4 = - 1.0 / 30.0
  double precision, parameter :: b6 =   1.0 / 42.0
  double precision, parameter :: b8 = - 1.0 / 30.0
!
  double precision trigamma
  double precision x
  double precision y
  double precision z
!
!  1): If X is not positive, fail.
!
  if ( x <= 0.0 ) then

    trigamma = 0.0
    write ( *, * ) ' '
    write ( *, * ) 'TRIGAMNA - Fatal error!'
    write ( *, * ) '  X <= 0.'
    stop
!
!  2): If X is smaller than A, use a small value approximation.
!
  else if ( x <= a ) then

    trigamma = 1.0 / x**2
!
!  3): Otherwise, increase the argument to ( X + I ) >= B...
!
  else

    z = x
    trigamma = 0.0

10      continue

    if ( z < b ) then
      trigamma = trigamma + 1.0 / z**2
      z = z + 1.0
      go to 10
    end if
!
!  ...and then apply an asymptotic formula.
!
    y = 1.0 / z**2

    trigamma = trigamma + 0.5 * &
           y + ( 1.0 &
         + y * ( b2 &
         + y * ( b4 &
         + y * ( b6 &
          + y *   b8 )))) / z

  end if

  return
end
function uniform_01_sample ( iseed )
!
!*******************************************************************************
!
!! UNIFORM_01_SAMPLE is a portable random number generator.
!
!
!  Formula:
!
!    ISEED = ISEED * (7**5) mod (2**31 - 1)
!    UNIFORM_01_SAMPLE = ISEED * / ( 2**31 - 1 )
!
!  Parameters:
!
!    Input/output, integer ISEED, the integer "seed" used to generate
!    the output random number, and updated in preparation for the
!    next one.  ISEED should not be zero.
!
!    Output, real UNIFORM_01_SAMPLE, a random value between 0 and 1.
!

!
!  IA = 7**5
!  IB = 2**15
!  IB16 = 2**16
!  IP = 2**31-1
!
  integer, parameter :: ia = 16807
  integer, parameter :: ib15 = 32768
  integer, parameter :: ib16 = 65536
  integer, parameter :: ip = 2147483647
!
  integer iprhi
  integer iseed
  integer ixhi
  integer k
  integer leftlo
  integer loxa
  double precision uniform_01_sample
!
!  Don't let ISEED be 0.
!
  if ( iseed == 0 ) then
    iseed = ip
  end if
!
!  Get the 15 high order bits of ISEED.
!
  ixhi = iseed / ib16
!
!  Get the 16 low bits of ISEED and form the low product.
!
  loxa = ( iseed - ixhi * ib16 ) * ia
!
!  Get the 15 high order bits of the low product.
!
  leftlo = loxa / ib16
!
!  Form the 31 highest bits of the full product.
!
  iprhi = ixhi * ia + leftlo
!
!  Get overflow past the 31st bit of full product.
!
  k = iprhi / ib15
!
!  Assemble all the parts and presubtract IP.  The parentheses are
!  essential.
!
  iseed = ( ( ( loxa - leftlo * ib16 ) - ip ) + ( iprhi - k * ib15 ) * ib16 ) &
    + k
!
!  Add IP back in if necessary.
!
  if ( iseed < 0 ) then
    iseed = iseed + ip
  end if
!
!  Multiply by 1 / (2**31-1).
!
  uniform_01_sample = dble ( iseed ) * 4.656612875e-10

  return
end
    
subroutine poisson(u,nc) !samples number of poissons from a poisson distribution 
!==================== 
INTEGER :: nc 
double precision :: u,x,y 
x=EXP(-1.0d0) !variance of a poission is equal to its mean, mean is 1 recombination per Morgan 
y=x 
nc=0 
DO 
IF(x>u) exit 
nc=nc+1 
y=y/dble(nc) 
x=x+y 
END do 
    END subroutine poisson 

 subroutine F_and_f_molecularAllo(nanim,nSNP,n_chr,geno,pOBS,fdev)

!     Computes molecular inbreeding from proportion of homozygous SNP (fmol1)
!     Computes coancestries (including self-coancestries) from identity of alleles
!     (as in Nejati-Javaremi et al. 1997)
!     Computes molecular inbreeding from self-coancestry as Fi=2fii-1 (fmol2)

!  Notes:
!  1. Includes computation of coancestry derived from expected homozygosity
!  2. Allows to correct by the frequency of parents
      

!  Variables:
!  coanc(i,j)= coancestry coefficient between animals i and j
!  count_miss1, count_miss2= counter for missing genotypes
!  fmol1= molecular F coefficient computed as proportion of homozygous SNP
!  fmol2= molecular F coefficient computed from self-coancestries
!  miss(i,j)= 9999 if genotype for SNP j of animal i is missing (0 otherwise)

!  B. Villanueva February, 2012
 
!!!=========================================================================================
!!! 
!!!Cambios para hacerlo subroutine
!!! Lo he cambiado para que en lugar de calcular las frecuencias de los padres lea el archivo de 
!!!    frecuencias de la generación 0 estimado en el programa. 
!!!Cambios de formato para coincidir con el programa principal

!!  ft=matriz de parentesco
!!  genot=genotipos
!!  kkkk= similutud(1) o L&H(2)
!!!
!!! Eli Enero, 2019
!!!=========================================================================================
 
!      use param1 !!!!Eli

implicit none
integer, intent(in) :: nanim,n_chr
integer, dimension(n_chr), intent(in) :: nSNP
integer*1, dimension (2,maxval(nSNP),n_chr,nanim),  intent(in)  :: geno
double precision :: z,mincoanc,maxcoanc,meancoanc
integer :: i,j,k,kk,l,icount,x,y,count_miss1,count_miss2,cont,nSNPt
double precision, dimension(sum(nSNP) ), intent(in) :: pOBS
integer, dimension (nanim) :: alelos
character(len=11) :: mark    
double precision, dimension(nanim,nanim),  intent(out) :: fdev    
integer, dimension(nanim) :: anim, par,nhomo    
integer, dimension(nanim,nanim) ::icount2
double precision, dimension(nanim) :: fmol1,fdev2,EHind    
double precision, dimension(nanim,nanim) :: coanc,EHtot,sumcoan    

integer*1, allocatable, dimension (:,:) :: genot     
integer, allocatable, dimension(:,:,:) :: allele
integer, allocatable, dimension(:,:) :: miss
double precision, allocatable, dimension(:) :: EH

    nSNPt=sum(nSNP)    
    
    
allocate(genot(nSNPt,nanim))   
allocate(allele(2,nSNPt,nanim))    
allocate(miss(nanim,nSNPt))    
!allocate(pOBS(nSNPt))
allocate(EH(nSNPt))     


     
genot=0   
kk=0
do j=1,n_chr
    !print*,j
    do l=1,nSNP(j)
        alelos=0
        kk=kk+1
        do i=1,nanim
            select case(sum(geno(:,l,j,i)))
            case(2)
                alelos(i)=2
            case(1)
                alelos(i)=1
            case(0)
                alelos(i)=0
            end select
        enddo
          
        genot(kk,:)=alelos(:) 

    enddo
enddo 

      

 
      

!     MOLECULAR
!     =========
      nhomo=0; allele=0
      coanc=0.0; fmol1=0.0
      miss=0; EHind=0.0; sumcoan=0.0
      fdev2=0.0; fdev=0.0 
          
      EH=0.0; EHtot=0.0; icount2=0



!open(43,file='SNP_freq_base.txt')
!
!      
! read(43,*)     
!
do k=1, nSNPt
     !read(43,*) pOBS(k)
    EH(k)=1-(2*pOBS(k)*(1-pOBS(k)))
    
end do    
!      
!close(43)
!      


        
! Calculating inbreeding coefficients    
   
      do i=1,nanim
        count_miss1=0

        do j=1,nSNPt

          if(genot(j,i)==0) then
            nhomo(i)=nhomo(i)+1
            allele(1,j,i)=1; allele(2,j,i)=1
          end if

          if(genot(j,i)==1) then
            allele(1,j,i)=1; allele(2,j,i)=2
          end if

          if(genot(j,i)==2) then
            nhomo(i)=nhomo(i)+1
            allele(1,j,i)=2; allele(2,j,i)=2
          end if

          if(genot(j,i)==3) then
            miss(i,j)=9999
            count_miss1=count_miss1+1
            allele(1,j,i)=0; allele(2,j,i)=0
          end if
          
          if (miss(i,j)/=9999) then
              EHind(i)=EHind(i)+EH(j)
          end if
        end do
        


        fmol1(i)=dble(nhomo(i))/dble(nSNPt-count_miss1)
        if(real(nSNPt-count_miss1-EHind(i)).eq.0) then 
        !print*, dble(nSNPt-count_miss1-EHind(i))
        endif
        fdev2(i)=dble(nhomo(i)-EHind(i))/dble(nSNPt-count_miss1-EHind(i))
        
      end do

      do i=1,nanim
        do j=1,nanim
          icount=0
          count_miss2=0
          do k=1,nSNPt
            if(miss(i,k)/=9999.and.miss(j,k)/=9999) then
              if(allele(1,k,i)==allele(1,k,j)) icount=icount+1
              if(allele(1,k,i)==allele(2,k,j)) icount=icount+1
              if(allele(2,k,i)==allele(1,k,j)) icount=icount+1
              if(allele(2,k,i)==allele(2,k,j)) icount=icount+1
              icount2(i,j)=icount2(i,j)+icount
              EHtot(i,j)=EHtot(i,j)+EH(k)
              else
              count_miss2=count_miss2+1
              end if
          end do
          

          
            coanc(i,j)=dble(icount)/dble(4*(nSNPt-count_miss2))
    !         PRINT*,i,j,anim(i),anim(j),coanc(i,j)
            sumcoan(i,j)=dble(icount)/4.00
            fdev(i,j)=dble(sumcoan(i,j)-EHtot(i,j))/dble((nSNPt-count_miss2-EHtot(i,j)))
            !PRINT*,i,j,EHtot(i,j),fdev(i,j)


          
        end do
      end do
      
      PRINT*, "Calculation of molecular coancestries finished"
      PRINT*, "ALL CALCULATIONS FINISHED AND SAVED"



      !deallocate(par)
      !deallocate(allele)
      !deallocate(nhomo)
      !
      !deallocate(miss)
      !deallocate(fmol1)
      !deallocate(sumcoan)
      !deallocate(pOBS)
      !deallocate(EH)
      !deallocate(coanc)
      !deallocate(EHtot)
      !deallocate(icount2)
      !!deallocate(fdev)
      !deallocate(fdev2)
      !deallocate(EHind)
      !deallocate(genot)  
      !
      


    end subroutine F_and_f_molecularAllo


 subroutine F_and_f_molecularAllo2(nanim,npob,ntotal,nSNP,n_chr,geno,pOBS,fdev)

!     Computes molecular inbreeding from proportion of homozygous SNP (fmol1)
!     Computes coancestries (including self-coancestries) from identity of alleles
!     (as in Nejati-Javaremi et al. 1997)
!     Computes molecular inbreeding from self-coancestry as Fi=2fii-1 (fmol2)

!  Notes:
!  1. Includes computation of coancestry derived from expected homozygosity
!  2. Allows to correct by the frequency of parents
      

!  Variables:
!  coanc(i,j)= coancestry coefficient between animals i and j
!  count_miss1, count_miss2= counter for missing genotypes
!  fmol1= molecular F coefficient computed as proportion of homozygous SNP
!  fmol2= molecular F coefficient computed from self-coancestries
!  miss(i,j)= 9999 if genotype for SNP j of animal i is missing (0 otherwise)

!  B. Villanueva February, 2012
 
!!!=========================================================================================
!!! 
!!!Cambios para hacerlo subroutine
!!! Lo he cambiado para que en lugar de calcular las frecuencias de los padres lea el archivo de 
!!!    frecuencias de la generación 0 estimado en el programa. 
!!!Cambios de formato para coincidir con el programa principal

!!  ft=matriz de parentesco
!!  genot=genotipos
!!  kkkk= similutud(1) o L&H(2)
!!!
!!! Eli Enero, 2019
!!!=========================================================================================
 
!      use param1 !!!!Eli

implicit none
integer, intent(in) :: nanim,n_chr,npob
integer, dimension(n_chr), intent(in) :: nSNP
integer*1, dimension (2,maxval(nSNP),n_chr,nanim),  intent(in)  :: geno
integer, dimension(npob,2), intent(in) :: ntotal
double precision :: z,mincoanc,maxcoanc,meancoanc,tt
integer :: i,j,k,t,ti,tj,o,kk,l,icount,x,y,count_miss1,count_miss2,cont,nSNPt
double precision, dimension(sum(nSNP),npob,npob), intent(in) :: pOBS
integer, dimension (nanim) :: alelos

character(len=11) :: mark    
double precision, dimension(nanim,nanim),  intent(out) :: fdev    
integer, dimension(nanim) :: anim, par,nhomo    
integer, dimension(nanim,nanim) ::icount2
double precision, dimension(nanim) :: fmol1,fdev2,EHind    
double precision, dimension(nanim,nanim) :: coanc,EHtot,sumcoan    

integer*1, allocatable, dimension (:,:) :: genot     
integer, allocatable, dimension(:,:,:) :: allele
integer, allocatable, dimension(:,:) :: miss
double precision, allocatable, dimension(:,:,:) :: EH
integer, allocatable, dimension(:) ::nindpop

    nSNPt=sum(nSNP)    
    

allocate(genot(nSNPt,nanim))   
allocate(allele(2,nSNPt,nanim))    
allocate(miss(nanim,nSNPt))    
!allocate(pOBS(nSNPt))
allocate(EH(nSNPt,npob,npob))
allocate(nindpop(npob)) 

     
genot=0   
kk=0
do j=1,n_chr
    !print*,j
    do l=1,nSNP(j)
        alelos=0
        kk=kk+1
        do i=1,nanim
            select case(sum(geno(:,l,j,i)))
            case(2)
                alelos(i)=2
            case(1)
                alelos(i)=1
            case(0)
                alelos(i)=0
            end select
        enddo
          
        genot(kk,:)=alelos(:) 

    enddo
enddo 

      

 
      

!     MOLECULAR
!     =========
      nhomo=0; allele=0
      coanc=0.0; fmol1=0.0
      miss=0; EHind=0.0; sumcoan=0.0
      fdev2=0.0; fdev=0.0 
          
      EH=0.0; EHtot=0.0; icount2=0



!open(43,file='SNP_freq_base.txt')
!
!      
! read(43,*)     
!
do i=1,npob
    do j=1,npob
        do k=1, nSNPt
     !read(43,*) pOBS(k)
            EH(k,i,j)=1-(2*pOBS(k,i,j)*(1-pOBS(k,i,j)))
        enddo
    end do    
end do
!      
!close(43)
!      


        
! Calculating inbreeding coefficients    
kk=1
o=0
i=1

   do k=1,npob
       o=o+(sum(ntotal(k,:)))
      do i=kk,o
        count_miss1=0

        do j=1,nSNPt

          if(genot(j,i)==0) then
            nhomo(i)=nhomo(i)+1
            allele(1,j,i)=1; allele(2,j,i)=1
          end if

          if(genot(j,i)==1) then
            allele(1,j,i)=1; allele(2,j,i)=2
          end if

          if(genot(j,i)==2) then
            nhomo(i)=nhomo(i)+1
            allele(1,j,i)=2; allele(2,j,i)=2
          end if

          if(genot(j,i)==3) then
            miss(i,j)=9999
            count_miss1=count_miss1+1
            allele(1,j,i)=0; allele(2,j,i)=0
          end if
          
          if (miss(i,j)/=9999) then
              EHind(i)=EHind(i)+EH(j,k,k)
          end if
        end do
        


        fmol1(i)=dble(nhomo(i))/dble(nSNPt-count_miss1)
        if(real(nSNPt-count_miss1-EHind(i)).eq.0) then 
        !print*, dble(nSNPt-count_miss1-EHind(i))
        endif
        fdev2(i)=dble(nhomo(i)-EHind(i))/dble(nSNPt-count_miss1-EHind(i))
    
      end do
          kk=kk+sum(ntotal(k,:))    
   end do

   
   kk=0
    do i=1,npob
        nindpop(i)=kk+sum(ntotal(i,:))
        kk=nindpop(i)
    enddo

      do i=1,nanim
        do j=1,nanim
            ti=0
            tj=0
            do t=1,npob
                if(t.eq.1.and.i.le.nindpop(t)) then
                    ti=t
                elseif(i.le.nindpop(t).and.i.gt.nindpop(t-1)) then
                    ti=t
                endif
             enddo
            do t=1,npob
                if(t.eq.1.and.j.le.nindpop(t)) then
                    tj=t
                elseif(j.le.nindpop(t).and.j.gt.nindpop(t-1)) then
                    tj=t
                endif
             enddo            
            
          icount=0
          count_miss2=0
          do k=1,nSNPt
            if(miss(i,k)/=9999.and.miss(j,k)/=9999) then
                if(allele(1,k,i)==allele(1,k,j)) icount=icount+1
                if(allele(1,k,i)==allele(2,k,j)) icount=icount+1
                if(allele(2,k,i)==allele(1,k,j)) icount=icount+1
                if(allele(2,k,i)==allele(2,k,j)) icount=icount+1
                icount2(i,j)=icount2(i,j)+icount
                tt=EH(k,ti,tj)                   
                EHtot(i,j)=EHtot(i,j)+tt
            else
              count_miss2=count_miss2+1
            end if
          end do
          

          
            coanc(i,j)=dble(icount)/dble(4*(nSNPt-count_miss2))
    !         PRINT*,i,j,anim(i),anim(j),coanc(i,j)
            sumcoan(i,j)=dble(icount)/4.00
            fdev(i,j)=dble(sumcoan(i,j)-EHtot(i,j))/dble((nSNPt-count_miss2-EHtot(i,j)))
            !PRINT*,i,j,EHtot(i,j),fdev(i,j)


          
        end do
      end do
      
      !PRINT*, "Calculation of molecular coancestries finished"
      !PRINT*, "ALL CALCULATIONS FINISHED AND SAVED"



      !deallocate(par)
      !deallocate(allele)
      !deallocate(nhomo)
      !
      !deallocate(miss)
      !deallocate(fmol1)
      !deallocate(sumcoan)
      !deallocate(pOBS)
      !deallocate(EH)
      !deallocate(coanc)
      !deallocate(EHtot)
      !deallocate(icount2)
      !!deallocate(fdev)
      !deallocate(fdev2)
      !deallocate(EHind)
      !deallocate(genot)  
      !
      


    end subroutine F_and_f_molecularAllo2

    


    
    
        
 subroutine Gmatrix(method,n_snp,n_chr,n_animal,genotype,freq,ft)

 !Archivo con genotipos codificados con 0, 1, 2 y 5 (NA)
 !Archivo con las frecuencias de los genotipos codificados con 0 en el archivo anterior

 ! modificaciones con respecto al original
 !  linea 91:
 !  do i =0,n_snp a do i =1,n_snp
 !  lineas
 ! se quitan las lineas que modifican la G para invertirla: lineas 145 a 148
 !  G=G*0.95
 !  do i=1,n_animal
 !     G(i,i)=G(i,i) + 0.05
 !  enddo
 ! modificado para que lo haga con todos los SNsPs no con los que tienen freq >0.05 modificada linea 104 y 131 >0 
 
 ! modificaciones ELI 
 !  añadi ft=matriz de parentesco
 !  añadi method 3=yang
 !  
 ! 
 !  
 ! modificado
 


 !cambia el orden de los datos en el archivo de salida

  implicit none
    !integer:: n_snp
    integer,parameter:: dp=kind(1d0)
    integer, intent(in):: method,n_animal,n_chr
    integer, dimension(n_chr), intent(in) :: n_snp
    double precision, dimension(sum(n_snp) ), intent(in) :: freq
    integer :: i,j,k,l,io,pos,coun,cont,rank,s, ms, A,nSNPt     
    integer*1, dimension (2,maxval(n_snp),n_chr,n_animal),  intent(in)  :: genotype
    integer, allocatable::geno(:,:),alelovr(:)
    character (len=5), allocatable::id(:)
    double precision, allocatable::X(:,:),G(:,:),XD(:,:),XP(:,:),T(:)
    double precision:: t1,t2,H
    double precision ::AX,X2,p2
    double precision, dimension(n_animal,n_animal) :: ft
    ft=0
    nSNPt=sum(n_snp) 


allocate(geno(n_animal,nSNPt))

do
!  print *,' 1 - VanRaden firstG = ZZ''/sum(2pq)'
!  print *,' 2 - VanRaden secondG = Yang et al. = mean (Z_i Z_i'' /(2p_i q_i))'

  !print*,'method of VanRaden'
  !read  *,method
  
  if(method==1 .or. method==2 .or. method==3) exit
enddo

!n_animal=0;io=0
   !write (*,*) ' reading data'
     call cpu_time(t1)

allocate(alelovr(nSNPt))     
     
    do i=1,n_animal
        alelovr=0
        !print*,i
        cont=0
        do j=1,n_chr
            do l=1,n_snp(j)     
                cont=cont+1
                select case(sum(genotype(:,l,j,i)))
                case(2) 
                    alelovr(cont)=2
                case(1)
                    alelovr(cont)=1
                case(0)
                    alelovr(cont)=0
                end select
                enddo
        enddo          
        geno(i,:)=alelovr(:)
    enddo 
     
     
DO  !Lee el n?mero de lineas en el archivo de datos
 !read (10,*,iostat=io)
 IF(io.ne.0) EXIT
 !n_animal=n_animal+1

! if (n_animal.eq.1000) exit
END DO
!rewind(10)
allocate (X(n_animal,nSNPt),XP(n_animal,nSNPt),&
          XD(n_animal,nSNPt),G(n_animal,n_animal),id(n_animal),T(n_animal))

   !write (*,*) n_animal,' individuals'
   !write (*,*) nSNPt,' snps'

  DO i=1,n_animal
    !read (10,*,iostat=io) id(i), geno(i,1:n_snp)
    !if (mod(i,1).eq.0.d0) !write (*,*) i,geno(i,1:nSNPt)

  END DO
  io=0
 call cpu_time(t2)
  !print *,'data read, time',t2-t1
  call cpu_time(t1)

  !print *,'X setup'
  X=0d0
  do i=1,nSNPt
     !if(mod(i,500)==0) !print *,i,'ith SNP'
     pos=i ! position of the first allele
     where(geno(:,i)==0)
        X(:,pos)=0d0
     elsewhere(geno(:,i)==1)
        X(:,pos)=1d0
     elsewhere(geno(:,i)==2)
        X(:,pos)=2d0
     elsewhere(geno(:,i)==5)
        X(:,pos)=1d0
     end where
  enddo
  !Calculo de frecuencias alélicas
  !write (*,*) ' calculating allele frequencies'

  ! deviate X from average within-individual value
  !open(11,file="SNP_freq_base.txt", form='formatted', status='old')

!read(11,*)
  do i =1,nSNPt
     !if(mod(i,1)==0d0) !print *,i,'ith SNP'
     !freq(i)=sum(X(:,i))/(2.d0*n_animal)
     !read(11,*) freq(i)
     !freq(i)=1-freq(i)
     XP(:,i) = X(:,i)-2.d0*freq(i)
  enddo
!pause
  write (*,*) ' calculating G elements'
!close(11)
! polymorphic SNPs
  coun=count( (freq*(1.d0-freq)) > 0)
  !print *,'coun',coun

  print *,'Computing G'


  !=================================================
  ! compute X D X' where D is diag(freq*(1-freq))
  ! Ignacio Aguilar - like style
  ! This is PVR "second G"
  ! too slow!!
  !=================================================
  !do i=1,nanim
  !  do j=1,nanim
  !    do k=1,nsnp
  !      if( (freq(k)*(1d0-freq(k))) >0) then
  !	      G(i,j) = G(i,j) + X(i,k)*X(j,k)/(2d0*freq(k)*(1d0-freq(k)))
  !      endif
  !    enddo
  !  enddo
!enddo
  
  
  print*, 'method',method

if(method==2 .or. method==3) then

    !------ Compute  VanRaden's second G 
    XD=0d0
    do i=1,nSNPt
        if( (freq(i)*(1d0-freq(i))) > 0) then
            XD(:,i)=XP(:,i)/sqrt(2.d0*freq(i)*(1.d0-freq(i)))
             
        endif
    enddo
    G=matmul(XD,transpose(XD))
    G=G/dble(coun)
elseif(method==1) then

    !------ Compute Van Raden first G
    G=matmul(XP,transpose(XP))/ (2.d0*sum(freq*(1d0-freq)))
endif

    if(method==3) then

        !------ Compute   Yang G

        do k=1,n_animal
            G(k,k)=0
            do i=1,nSNPt
                if( (freq(i)*(1d0-freq(i))) > 0) then       
                AX=(1+(2*(freq(i))))*(X(k,i))
                X2=(X(k,i))*(X(k,i))
                p2=2*((freq(i))*(freq(i)))
                T(k)=T(k)+((X2-AX+p2)/(2*freq(i)*(1-freq(i))))
                endif

            enddo
                G(k,k)=1+(T(k)/coun)

        enddo
        
    endif
    
   
     

  call cpu_time(t2)
  !print *,'G computed, time',t2-t1


  
  ft=dble(G)



print*, 'DONE'

end subroutine Gmatrix
 
    
    
       
 subroutine Gmatrix3(method,npob,ntotal,n_snp,n_chr,n_animal,genotype,freq,ft)

 !Archivo con genotipos codificados con 0, 1, 2 y 5 (NA)


  implicit none
    !integer:: n_snp
    integer,parameter:: dp=kind(1d0)
    integer, intent(in):: method,n_animal,n_chr,npob
    integer, dimension(n_chr), intent(in) :: n_snp
    integer, dimension(npob,2), intent(in) :: ntotal
    double precision, dimension(sum(n_snp),npob,npob), intent(in) :: freq
    double precision, dimension(sum(n_snp),npob,npob):: EH
    integer :: i,j,jj,k,k2,kk,lj,lk,l,io,pos,coun,cont,rank,s, ms, A,nSNPt     
    integer*1, dimension (2,maxval(n_snp),n_chr,n_animal),  intent(in)  :: genotype
    integer, allocatable::geno(:,:),alelovr(:)
    double precision, allocatable::X(:,:),G(:,:),XP(:,:),T(:)
    double precision:: t1,t2,H,xx
    double precision ::AX,X2,p2
    double precision, dimension(n_animal,n_animal) :: ft
    integer, allocatable, dimension(:) ::nindpop
    ft=0
    nSNPt=sum(n_snp) 
    
allocate(geno(n_animal,nSNPt))
allocate(nindpop(npob)) 
geno=0
nindpop=0

  !print *,' 1 - VanRaden firstG = ZZ''/sum(2pq)'
  !print *,' 2 - VanRaden secondG = Yang et al. = mean (Z_i Z_i'' /(2p_i q_i))'

  print*,'method of VanRaden'
  !read  *,method
  


!n_animal=0;io=0
   !write (*,*) ' reading data'
     call cpu_time(t1)

allocate(alelovr(nSNPt))     
alelovr=0
     
    do i=1,n_animal
        alelovr=0
        !print*,i
        cont=0
        do j=1,n_chr
            do l=1,n_snp(j)     
                cont=cont+1
                select case(sum(genotype(:,l,j,i)))
                case(2) 
                    alelovr(cont)=2
                case(1)
                    alelovr(cont)=1
                case(0)
                    alelovr(cont)=0
                end select
                enddo
        enddo          
        geno(i,:)=alelovr(:)
    enddo 
 

!rewind(10)
allocate (X(n_animal,nSNPt),XP(n_animal,n_animal),&
          G(n_animal,n_animal),T(n_animal))
X=0
XP=0
G=0
T=0

   write (*,*) n_animal,' individuals'
   write (*,*) nSNPt,' snps'


  io=0
 
 call cpu_time(t2)
  print *,'data read, time',t2-t1
  call cpu_time(t1)

  print *,'X setup'
  X=0d0

  x=geno
  


!read(11,*)
  EH=0
  
  do i=1,npob
    do j=1,npob
        do k=1,nSNPt
            EH(k,i,j)=1-(2*freq(k,i,j)*(1-freq(k,i,j)))
        enddo
    end do    
end do

   kk=0
    do i=1,npob
        nindpop(i)=kk+sum(ntotal(i,:))
        kk=nindpop(i)
    enddo
    
    
    print*, 'method=',method
if(method==2 .or. method==3) then
    
    
    
do j=1,n_animal
   do k=1,n_animal
     do l=1,npob
        if(l.eq.1.and.j.le.nindpop(l)) then
            lj=l
        elseif(l.gt.1.and.j.le.nindpop(l).and.j.gt.nindpop(l-1)) then
            lj=l
        endif

        if(l.eq.1.and.k.le.nindpop(l)) then
            lk=l
        elseif(l.gt.1.and.k.le.nindpop(l).and.k.gt.nindpop(l-1)) then
            lk=l
        endif
    enddo      
    do i =1,nSNPt
     !if(mod(i,1)==0d0) print *,i,'ith SNP'
     xx = (X(j,i)-2.d0*freq(i,lj,lk))*(X(k,i)-2.d0*freq(i,lj,lk))/EH(i,lj,lk)
     XP(j,k) = XP(j,k)+xx
    enddo
    G(j,k)=(1/2*nSNPt)*(XP(j,k))
   enddo
enddo

endif


if(method==3) then

    !------ Compute   Yang G

    do k=1,n_animal
        G(k,k)=0
        do l=1,npob
            if(l.eq.1.and.k.le.nindpop(l)) then
                j=l
            elseif(k.le.nindpop(l).and.j.gt.nindpop(l-1)) then
                j=l
            endif
        enddo
        do i=1,nSNPt
            if( (freq(i,j,j)*(1d0-freq(i,j,j))) > 0) then       
            AX=(1+(2*(freq(i,j,j))))*(X(k,i))
            X2=(X(k,i))*(X(k,i))
            p2=2*((freq(i,j,j))*(freq(i,j,j)))
            T(k)=T(k)+((X2-AX+p2)/(2*freq(i,j,j)*(1-freq(i,j,j))))
            endif
        enddo
        G(k,k)=1+(T(k)/nSNPt)
    enddo
endif
    
  
  call cpu_time(t2)
  print *,'G computed, time',t2-t1
 
  ft=dble(G)

print*, 'DONE'

end subroutine Gmatrix3
 



  