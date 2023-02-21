module variables
implicit none
integer :: n_chr,genom_length,n_neut_total,n_neut,n_SNP,n_SNP_total,n_pos,n_posf,nparent,&
        SNP_pos_nofix,neut_pos_nofix,n_SNPf,n_neutf,SNP_pos,neut_pos,n_gen,aux,auxse,&
        SNP_posn,neut_posn,n_SNPn,n_neutn,n_posn,asx,asxse,atx,atxse,n_SNP_total_new,n_neut_total_new
integer :: n_ind,nbase,genom_map,genom_map_neut,genom_map_nofix,genom_map_neut_nofix,&
        n_genopt,i,j,k,l,m,n,n_SNP_fix,n_neut_fix,segrega,jcr,gen_1,gen_2,&
        genom_map_nuevo,genom_map_neut_nuevo

integer*1 genotype,genotype_nofix,genotype_t1,genotype_nuevo
integer*8 total_genom_length
double precision :: results,f_mol_neut,f_mol_SNP,frec,time1,time2,inb_mol,coanc_mol,&
        exp_het,obs_het,freq_nM,r2,dprime,ld,ld2
real :: u
character*3 arep


dimension :: n_ind(2),inb_mol(3,2),coanc_mol(3,2),&           
          exp_het(3,2),obs_het(3,2),segrega(2)
allocatable :: genom_length(:),n_neut(:),n_SNP(:),n_pos(:),n_posf(:),f_mol_SNP(:,:),frec(:,:),results(:,:),&
            genotype(:,:,:,:),genom_map(:,:),genom_map_neut(:,:),genom_map_nofix(:,:),genotype_t1(:,:,:,:),&
			genom_map_neut_nofix(:,:),SNP_pos(:,:),neut_pos(:,:),f_mol_neut(:,:),genotype_nofix(:,:,:,:),aux(:),auxse(:),&
            genom_map_nuevo(:,:),genom_map_neut_nuevo(:,:),genotype_nuevo(:,:,:,:),asx(:,:),asxse(:,:),n_posn(:),&
			atx(:,:),atxse(:,:)
allocatable :: gen_1(:),gen_2(:)
allocatable :: SNP_pos_nofix(:,:),neut_pos_nofix(:,:),n_SNPf(:),n_neutf(:)
allocatable :: SNP_posn(:,:),neut_posn(:,:),n_SNPn(:),n_neutn(:)

end module

program filtrobase
use variables
implicit none
integer :: s,t,iseed,cont,m2,op,m3,m4,m5, count_rate, count_max,matriz,diez,son  !!,optt,rl2,rl
integer :: igs,ngs,n_cross,cross,pat,sobant,geneal,oploci
character(len=11) :: a
double precision :: long,long0,interval
double precision :: maf,varmaf,mafneut,varmafneut,freq_filtro,freq_filtrob
double precision, allocatable, dimension(:) :: x, qx, xn, qxn
double precision, allocatable, dimension(:) :: pOBS,pOBSneut,EH,freqmaf,freqmafneut
!character, allocatable, dimension(:,:,:,:) :: gk

allocatable cross(:),geneal(:,:),interval(:)
allocate(cross(1))

!call CPU_TIME(iseed)
!CALL SYSTEM_CLOCK(iseed, count_rate, count_max)
open(UNIT=55,FILE='repl.txt')
read(55,*) iseed
close(55)

call setup(iseed)
open(UNIT=55,FILE='filtrorepl_inputv4.txt')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read from file the parameters for simulations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! No. individuals (males and females, no. individuals total, no. chromosomes !
read(55,*) n_ind,nbase,n_chr
       nparent=sum(n_ind)


print*,n_ind,nparent,n_chr

allocate(n_neut(n_chr),n_SNP(n_chr),n_SNPf(n_chr),n_pos(n_chr),n_posf(n_chr),genom_length(n_chr),results(0:1,51))

n_neut=0
n_SNP=0
n_SNPf=0
n_pos=0
n_posf=0
genom_length=0
results=0

! No. total of SNP and non-marker loci !
read(55,*) n_SNP_total,n_neut_total
print*,n_SNP_total,n_neut_total





! total length of genome !
    read(55,*) total_genom_length
    	print*,total_genom_length
        
        
!!!Choose the filtering frequency op=1 remove snp fixed, op=2 remove snp and no-markers with freq_filtro,op=3 remove snp with freqfiltro and no-markers fixed!
        
read(55,*) op
print*,op


!!!Choose the filtering frequency
read(55,*) freq_filtro

!!!!!!Filter number of loci
read(55,*) oploci
read(55,*) n_SNP_total_new,n_neut_total_new
print*,oploci,n_SNP_total_new,n_neut_total_new

!!!!!Write ped with all loci optt=2, or write ped only with SNPs optt=1

!read(55,*) optt
!read(55,*)matriz

!  N. generaciones sin control, 
read(55,*) ngs

close(55)



! Total no. positions in each chromosome (SNP + non markers) !      
! Read  a file the genome maps characteristics !


open(UNIT=35,FILE='../../../maps.txt',status='old')

read(35,*)
do i=1,n_chr
 read(35,*) a,jcr,genom_length(i),n_SNP(i),n_neut(i)
enddo
allocate(SNP_pos(n_chr,maxval(n_SNP)))
allocate(genom_map(n_chr,maxval(n_SNP)))

SNP_pos=0
genom_map=0

!print*,'maxval(n_SNP)=', maxval(n_SNP)
n_pos=n_neut+n_SNP
n_posf=n_pos
print*,n_pos


read(35,*)
do i=1,n_chr 
  read(35,*)
  do j=1,n_SNP(i)
    read(35,*) SNP_pos(i,j),genom_map(i,j)
  enddo
enddo

close(35)

open(UNIT=35,FILE='../../../mapsneut.txt',status='old')


allocate(neut_pos(n_chr,maxval(n_neut)))
allocate(genom_map_neut(n_chr,maxval(n_neut)))

neut_pos=0
genom_map_neut=0



do i=1,n_chr 
  read(35,*)
  do j=1,n_neut(i)
    read(35,*) neut_pos(i,j),genom_map_neut(i,j)
  enddo
enddo

close(35)



! Genotypes of the base population !
! Read the initial genotypes of the base population !

!allocate(genotype_base(2,maxval(n_pos),n_chr,nbase))
!genotype_base=0
allocate(aux(nbase),auxse(sum(n_ind)))
aux=0
auxse=0

! Sample individuals from the base population !
  do i=1,nbase
    aux(i)=i
  enddo
  
  call desordena(aux,nbase)
  
  
  
do i=1,(sum(n_ind))
    auxse(i)=aux(i)
enddo

call ordena(auxse,sum(n_ind))




  allocate(genotype(2,maxval(n_pos),n_chr,sum(n_ind)))
    genotype=0
  

open(UNIT=25,FILE='../../../expanbase.txt',status='old')
m=1
    do i=1,nbase
        print*,i       
        if(i.eq.auxse(m)) then
        do j=1,n_chr
            read(25,*) k,l,((genotype(k,l,j,m),k=1,2),l=1,n_pos(j))
        enddo
        m=min(m+1,sum(n_ind))
        else
            do j=1,n_chr
                read(25,*) 
            enddo
        endif
    enddo


    close(25)

deallocate(aux,auxse)
  
    
 ! Calculating observed  allele frequencies and filtered loci fixed     
  PRINT*, "Calculating allele frequencies"   
        
allocate(pOBS(sum(n_SNP)))
allocate(pOBSneut(sum(n_neut)))
pOBS=0
pOBSneut=0       

!open(11,file="SNP_freq_MAF0_segr.txt")
!open(111,file="SNP_freq_base_segr.txt")
!open(12,file="neut_freq_MAF0_segr.txt")
!open(112,file="neut_freq_base_segr.txt")

call frequency(nparent,n_chr,n_pos,n_SNP,n_neut,SNP_pos,genotype,pOBS,pOBSneut)
!print*, 'freq', pOBS
!print*, 'freqneut', pOBSneut
                 

freq_filtrob=1-freq_filtro


select case(op)
    case(1)
    
        l=0
    do j=1,n_chr
        n_SNP_fix=0
        do k=1,n_SNP(j)    
            l=l+1
            if(pOBS(l).le.freq_filtro.or.pOBS(l).ge.freq_filtrob) then
                !print*,j,l,k
                !print*,pOBS(l)
                n_SNP_fix=n_SNP_fix+1
            endif
        enddo
        print*,'n_SNP_fix=', n_SNP_fix
        n_SNPf(j)=(n_SNP(j))-n_SNP_fix
        n_posf(j)=(n_posf(j))-n_SNP_fix
    enddo 
    
        allocate(genotype_nofix(2,maxval(n_posf),n_chr,nparent))    
        allocate(SNP_pos_nofix(n_chr,maxval(n_SNPf)),freqmaf(sum(n_SNPf)),freqmafneut(sum(n_neut)))
        allocate(genom_map_nofix(n_chr,maxval(n_SNPf)))
        !l=maxval(n_SNPf)
        !print*, 'lmaxval=',l
        !print*, n_SNPf

        genotype_nofix=0
        SNP_pos_nofix=0

        l=0
        n=0
        maf=0
        mafneut=0
        

        do j=1,n_chr
            m2=0
            m=0        
            m3=0
            do k=1,n_pos(j)
                m3=m3+1
                if(k.ne.SNP_pos(j,m3)) then 
                    m=m+1
                    n=n+1
                    m3=m3-1
                    genotype_nofix(:,m,j,:)=genotype(:,k,j,:)
                    !write(112,*) pOBSneut(n)
                    if(pOBSneut(n).GT.0.5) then 
                        pOBSneut(n)=1-pOBSneut(n)
                        mafneut=mafneut+pOBSneut(n)
                    end if
                    !write(12,*) pOBSneut(n)
                    freqmafneut(n)=pOBSneut(n)
                else      
                    l=l+1
                    if(pOBS(l).gt.freq_filtro.and.pOBS(l).lt.freq_filtrob) then 
                        m=m+1 
                        m2=m2+1
                        SNP_pos_nofix(j,m2)=m       
                        genotype_nofix(:,m,j,:)=genotype(:,k,j,:)          
                        genom_map_nofix(j,m2)=genom_map(j,m3)
                       ! write(111,*) pOBS(l)
                        if(pOBS(l).GT.0.5) then 
                            pOBS(l)=1-pOBS(l)
                            maf=maf+pOBS(l)
                        end if
                       ! write(11,*) pOBS(l)
                        freqmaf(m2)=pOBS(l)
                    endif
                endif
            enddo 
        end do    

    case(2)
        
        l=0
        do j=1,n_chr
            n_SNP_fix=0
            do k=1,n_SNP(j)    
                l=l+1
                if(pOBS(l).le.freq_filtro.or.pOBS(l).ge.freq_filtrob) then
                    !print*,j,l,k
                    !print*,pOBS(l)
                    n_SNP_fix=n_SNP_fix+1
                endif
            enddo
            print*,'n_SNP_fix=', n_SNP_fix
            n_SNPf(j)=(n_SNP(j))-n_SNP_fix
            n_posf(j)=(n_posf(j))-n_SNP_fix
        enddo 
        
        allocate(n_neutf(n_chr))
        n_neutf=0
        l=0
        
        do j=1,n_chr
            n_neut_fix=0
            do k=1,n_neut(j)    
                l=l+1
                if(pOBSneut(l).le.freq_filtro.or.pOBSneut(l).ge.freq_filtrob) then
                   
                    n_neut_fix=n_neut_fix+1
                endif
            enddo
            n_neutf(j)=(n_neut(j))-n_neut_fix
            n_posf(j)=(n_posf(j))-n_neut_fix
        enddo 
        
        allocate(genotype_nofix(2,maxval(n_posf),n_chr,nparent))    
        allocate(SNP_pos_nofix(n_chr,maxval(n_SNPf)),freqmaf(sum(n_SNPf)),freqmafneut(sum(n_neutf)))
        allocate(neut_pos_nofix(n_chr,maxval(n_neutf)))
        allocate(genom_map_nofix(n_chr,maxval(n_SNPf)))
        allocate(genom_map_neut_nofix(n_chr,maxval(n_neutf)))


        genotype_nofix=0
        SNP_pos_nofix=0
        neut_pos_nofix=0
        genom_map_neut_nofix=0
        genom_map_nofix=0
        l=0
        n=0
        maf=0
        mafneut=0


        do j=1,n_chr
            m2=0
            m3=1
            m=0
            m4=0
            m5=0
            do k=1,n_pos(j)

                !if(m3.le.n_SNP(j)) then
                    if(k.ne.SNP_pos(j,m3)) then 
                        n=n+1
                        m5=min(m5+1,n_neut(j))
                        if(pOBSneut(n).gt.freq_filtro.and.pOBSneut(n).lt.freq_filtrob) then 
                            m=m+1
                            m4=m4+1
                            genotype_nofix(:,m,j,:)=genotype(:,k,j,:)
                            neut_pos_nofix(j,m4)=m
                            genom_map_neut_nofix(j,m4)=genom_map_neut(j,m5)
                            !write(112,*) pOBSneut(n)
                            
                            if(pOBSneut(n).GT.0.5) then 
                                pOBSneut(n)=1-pOBSneut(n)
                                mafneut=mafneut+pOBSneut(n)
                            end if
                            !write(12,*) pOBSneut(n)
                            freqmafneut(m4)=pOBSneut(n)
                        end if
                    else      
                        l=l+1
                        if(pOBS(l).gt.freq_filtro.and.pOBS(l).lt.freq_filtrob) then 
                            m=m+1 
                            m2=m2+1
                            SNP_pos_nofix(j,m2)=m       
                            genotype_nofix(:,m,j,:)=genotype(:,k,j,:)          
                            genom_map_nofix(j,m2)=genom_map(j,m3)
                            !write(111,*) pOBS(l)
                            
                            if(pOBS(l).GT.0.5) then 
                                pOBS(l)=1-pOBS(l)
                                maf=maf+pOBS(l)
                            end if
                            !write(11,*) pOBS(l)
                            freqmaf(m2)=pOBS(l)
                        endif
                        m3=min(m3+1,n_SNP(j))
                    endif

            enddo 
        end do    
        
        print*,'n_neut', n_neut
        n_neut=n_neutf
        n_neut_total=sum(n_neut)
        deallocate(n_neutf)
        deallocate(neut_pos)
        deallocate(genom_map_neut)

        allocate(neut_pos(n_chr,maxval(n_neut)))
        allocate(genom_map_neut(n_chr,maxval(n_neut)))
        neut_pos=0
        genom_map_neut=0
        neut_pos=neut_pos_nofix
        genom_map_neut=genom_map_neut_nofix

        deallocate(neut_pos_nofix)
        deallocate(genom_map_neut_nofix)      
    case(3)
        
        l=0
        do j=1,n_chr
            n_SNP_fix=0
            do k=1,n_SNP(j)    
                l=l+1
                if(pOBS(l).le.freq_filtro.or.pOBS(l).ge.freq_filtrob) then
                    !print*,j,l,k
                    !print*,pOBS(l)
                    n_SNP_fix=n_SNP_fix+1
                endif
            enddo
            print*,'n_SNP_fix=', n_SNP_fix
            n_SNPf(j)=(n_SNP(j))-n_SNP_fix
            n_posf(j)=(n_posf(j))-n_SNP_fix
        enddo 
        
        allocate(n_neutf(n_chr))
        n_neutf=0
        l=0
        
        do j=1,n_chr
            n_neut_fix=0
            do k=1,n_neut(j)    
                l=l+1
                if(pOBSneut(l).le.0.or.pOBSneut(l).ge.1) then
                   
                    n_neut_fix=n_neut_fix+1
                endif
            enddo
            n_neutf(j)=(n_neut(j))-n_neut_fix
            n_posf(j)=(n_posf(j))-n_neut_fix
        enddo 
        
        allocate(genotype_nofix(2,maxval(n_posf),n_chr,nparent))    
        allocate(SNP_pos_nofix(n_chr,maxval(n_SNPf)),freqmaf(sum(n_SNPf)),freqmafneut(sum(n_neutf)))
        allocate(neut_pos_nofix(n_chr,maxval(n_neutf)))
        allocate(genom_map_nofix(n_chr,maxval(n_SNPf)))
        allocate(genom_map_neut_nofix(n_chr,maxval(n_neutf)))


        genotype_nofix=0
        SNP_pos_nofix=0
        neut_pos_nofix=0
        genom_map_neut_nofix=0
        genom_map_nofix=0
        l=0
        n=0
        maf=0
        mafneut=0

        do j=1,n_chr
            m2=0
            m3=1
            m=0
            m4=0
            m5=0
            do k=1,n_pos(j)

                !if(m3.le.n_SNP(j)) then
                    if(k.ne.SNP_pos(j,m3)) then 
                        n=n+1
                        m5=min(m5+1,n_neut(j))
                        if(pOBSneut(n).gt.0.and.pOBSneut(n).lt.1) then 
                            m=m+1
                            m4=m4+1
                            genotype_nofix(:,m,j,:)=genotype(:,k,j,:)
                            neut_pos_nofix(j,m4)=m
                            genom_map_neut_nofix(j,m4)=genom_map_neut(j,m5)
                            !write(112,*) pOBSneut(n)
                            
                            if(pOBSneut(n).GT.0.5) then 
                                pOBSneut(n)=1-pOBSneut(n)
                                mafneut=mafneut+pOBSneut(n)
                            end if
                            !write(12,*) pOBSneut(n)
                            freqmafneut(m4)=pOBSneut(n)
                        end if
                    else      
                        l=l+1
                        if(pOBS(l).gt.freq_filtro.and.pOBS(l).lt.freq_filtrob) then 
                            m=m+1 
                            m2=m2+1
                            SNP_pos_nofix(j,m2)=m       
                            genotype_nofix(:,m,j,:)=genotype(:,k,j,:)          
                            genom_map_nofix(j,m2)=genom_map(j,m3)
                            !write(111,*) pOBS(l)
                            
                            if(pOBS(l).GT.0.5) then 
                                pOBS(l)=1-pOBS(l)
                                maf=maf+pOBS(l)
                            end if
                            !write(11,*) pOBS(l)
                            freqmaf(m2)=pOBS(l)
                        endif
                        m3=min(m3+1,n_SNP(j))
                    endif
            enddo 
        end do    
        
        print*,'n_neut', n_neut
        n_neut=n_neutf
        n_neut_total=sum(n_neut)
        deallocate(n_neutf)
        deallocate(neut_pos)
        deallocate(genom_map_neut)

        allocate(neut_pos(n_chr,maxval(n_neut)))
        allocate(genom_map_neut(n_chr,maxval(n_neut)))
        neut_pos=0
        genom_map_neut=0
        neut_pos=neut_pos_nofix
        genom_map_neut=genom_map_neut_nofix

        deallocate(neut_pos_nofix)
        deallocate(genom_map_neut_nofix)             
        
end select
    





n_SNP=n_SNPf
n_pos=n_posf
n_SNP_total=sum(n_SNP)


deallocate(n_SNPf)
deallocate(genotype)
deallocate(SNP_pos)
deallocate(genom_map)

allocate(genotype(2,maxval(n_pos),n_chr,nparent))
allocate(SNP_pos(n_chr,maxval(n_SNP)))
allocate(genom_map(n_chr,maxval(n_SNP)))
SNP_pos=0

genotype=0
genom_map=0

SNP_pos=SNP_pos_nofix
genotype=genotype_nofix
genom_map=genom_map_nofix
             

deallocate(genotype_nofix)
deallocate(SNP_pos_nofix)
deallocate(genom_map_nofix)

!close(11)
!close(111) 
!close(12)
!close(112) 





!!Filtro para numero de loci
select case(oploci)
    case(1)
    
    
        deallocate(pOBS)
        deallocate(pOBSneut)
        allocate(pOBS(sum(n_SNP)))
        allocate(pOBSneut(sum(n_neut)))
        pOBS=0
        pOBSneut=0       
        call frequency(nparent,n_chr,n_pos,n_SNP,n_neut,SNP_pos,genotype,pOBS,pOBSneut)

        allocate(n_posn(n_chr),n_SNPn(n_chr),n_neutn(n_chr))

        n_SNPn=int(n_SNP_total_new/n_chr)
        n_neutn=int(n_neut_total_new/n_chr)
        n_posn=n_SNPn+n_neutn
        allocate(asx(n_chr,maxval(n_SNP)),asxse(n_chr,maxval(n_SNPn)))
        allocate(atx(n_chr,maxval(n_neut)),atxse(n_chr,maxval(n_neutn)))
        
        asx=0
        asxse=0       
        
        do j=1,n_chr
            m3=1
            m4=0
            do i=1,n_pos(j)
                if(i.ne.SNP_pos(j,m3)) then
                    m4=m4+1
                    atx(j,m4)=i
                else
                    asx(j,m3)=i
                    m3=min(m3+1,n_SNP(j))
                endif
            enddo
            call desordena(asx(j,:),n_SNP(j))
            call desordena(atx(j,:),n_neut(j))   
        enddo
          
  
        do j=1,n_chr
            m3=1
            m4=0
            m5=0
            do i=1,max(n_SNPn(j),n_neutn(j)) 
                if(n_SNPn(j).eq.n_neutn(j)) then
                    asxse(j,i)=asx(j,i)
                    atxse(j,i)=atx(j,i)
                elseif(n_SNPn(j).ne.n_neutn(j)) then
                    if(i.le.n_SNPn(j).and.i.le.n_neutn(j)) then
                        asxse(j,i)=asx(j,i)
                        atxse(j,i)=atx(j,i)
                    elseif(i.gt.n_SNPn(j).and.i.le.n_neutn(j))then
                        atxse(j,i)=atx(j,i)
                    elseif(i.gt.n_neutn(j).and.i.le.n_SNPn(j))then
                        asxse(j,i)=asx(j,i)
                    endif
                endif
            enddo
            call ordena(asxse(j,:),n_SNPn(j))
            call ordena(atxse(j,:),n_neutn(j))
        enddo

        deallocate(freqmaf,freqmafneut)
        
        allocate(genotype_nuevo(2,maxval(n_posn),n_chr,nparent))    
        allocate(SNP_posn(n_chr,maxval(n_SNPn)),freqmaf(sum(n_SNPn)),freqmafneut(sum(n_neutn)))
        allocate(neut_posn(n_chr,maxval(n_neutn)))
        allocate(genom_map_nuevo(n_chr,maxval(n_SNPn)))
        allocate(genom_map_neut_nuevo(n_chr,maxval(n_neutn)))
        genotype_nuevo=0
        SNP_posn=0
        neut_posn=0
        genom_map_neut_nuevo=0
        genom_map_nuevo=0
        l=0
        n=0
        s=0
        t=0
        maf=0
        mafneut=0
        freqmaf=0
        freqmafneut=0

        do j=1,n_chr
            m2=1
            m3=1
            m=0
            m4=1
            m5=0
            do k=1,n_pos(j)

                !if(m3.le.n_SNP(j)) then
                    if(k.ne.SNP_pos(j,m3)) then 
                        n=n+1
                        !m5=min(m5+1,n_neut(j))
                        m5=m5+1
                        if((k.eq.atxse(j,m4))) then 
                            m=m+1
                            s=s+1
                            genotype_nuevo(:,m,j,:)=genotype(:,k,j,:)
                            neut_posn(j,m4)=m
                            genom_map_neut_nuevo(j,m4)=genom_map_neut(j,m5)
                            !write(112,*) pOBSneut(n)
                            
                            if(pOBSneut(n).GT.0.5) then 
                                pOBSneut(n)=1-pOBSneut(n)
                                mafneut=mafneut+pOBSneut(n)
                            end if
                            !write(12,*) pOBSneut(n)
                            freqmafneut(s)=pOBSneut(n)
                            m4=min(m4+1,n_neutn(j))
                        end if
                    else      
                        l=l+1
                        if((k.eq.asxse(j,m2))) then  
                            m=m+1
                            t=t+1
                            SNP_posn(j,m2)=m       
                            genotype_nuevo(:,m,j,:)=genotype(:,k,j,:)          
                            genom_map_nuevo(j,m2)=genom_map(j,m3)
                            !write(111,*) pOBS(l)
                            
                            if(pOBS(l).GT.0.5) then 
                                pOBS(l)=1-pOBS(l)
                                maf=maf+pOBS(l)
                            end if
                            !write(11,*) pOBS(l)
                            freqmaf(t)=pOBS(l)
                            m2=min(m2+1,n_SNPn(j))
                        endif
                        m3=min(m3+1,n_SNP(j))
                        !m3=m3+1
                    endif

            enddo 
        end do    
        
        
        deallocate(asx,asxse)
        n_SNP=n_SNPn
        n_SNP_total=sum(n_SNP)
        deallocate(n_SNPn)
        deallocate(SNP_pos)
        deallocate(genom_map)
        allocate(SNP_pos(n_chr,maxval(n_SNP)))
        allocate(genom_map(n_chr,maxval(n_SNP)))
        SNP_pos=0
        genom_map=0
        SNP_pos=SNP_posn
        genom_map=genom_map_nuevo
        deallocate(SNP_posn)
        deallocate(genom_map_nuevo)  
        
        n_neut=n_neutn
        n_neut_total=sum(n_neut)
        deallocate(n_neutn)
        deallocate(neut_pos)
        deallocate(genom_map_neut)
        allocate(neut_pos(n_chr,maxval(n_neut)))
        allocate(genom_map_neut(n_chr,maxval(n_neut)))
        neut_pos=0
        genom_map_neut=0
        neut_pos=neut_posn
        genom_map_neut=genom_map_neut_nuevo
        deallocate(neut_posn)
        deallocate(genom_map_neut_nuevo)      
        
        n_pos=n_SNP+n_neut
        
        
end select

deallocate(pOBS)
deallocate(pOBSneut)

allocate(x(sum(n_SNP)))
allocate(qx(sum(n_SNP)))
maf=maf/sum(n_SNP)
varmaf=0
varmafneut=0

do i=1,sum(n_SNP)
    x(i)= freqmaf(i)-maf
    qx(i)=(x(i))*(x(i))
enddo
varmaf=(sum(qx))/sum(n_SNP)


allocate(xn(sum(n_neut)))
allocate(qxn(sum(n_neut)))
mafneut=mafneut/sum(n_neut)

do i=1,sum(n_neut)
    xn(i)= freqmafneut(i)-mafneut
    qxn(i)=(xn(i))*(xn(i))
enddo
varmafneut=(sum(qxn))/sum(n_neut)




            allocate(geneal(sum(n_ind),2),interval(n_chr))
            interval=dble(genom_length)/dble(n_pos)
    ! Bucle para generaciones sin control*
      do igs=1,ngs
          

            ! Decides (randomly) the parents of next generation !
              do i=1,sum(n_ind)
                do k=1,2
                  n=(k-1)*n_ind(1)
	              call mother(u)
	              l=min(int(u*n_ind(k))+1+n,n_ind(k)+n)
	              geneal(i,k)=l
	            enddo
              enddo

              allocate(genotype_t1(2,maxval(n_pos),n_chr,sum(n_ind)))

            ! Generates the new population !
              do j=1,n_chr
                long0=dble(genom_length(j))/100000000
                long=long0
                do i=1,sum(n_ind)
                  do k=1,2

            ! Determines the number and position !
            ! of crossovers and sort them        !

            ! For small chromosomes (less than 1 Morgan) !
                    if(long0.lt.1) then
		              call mother(u)
		              if(u.lt.long0) then
		                long=1
		              else
		                long=0
                      endif
                    else
                    long=long0
		            endif
		            n_cross=0
		            do n=1,anint(long)
		              call mother(u)
                      u=min(u,.999999)
		              call poisson(u,l) 
                      n_cross=n_cross+l
		            enddo
                    deallocate (cross)
                    allocate (cross(n_cross))
        
                    do n=1,n_cross
                      call mother(u)
                      cross(n)=min(int(u*genom_length(j)+1),genom_length(j))
		              cross(n)=anint(cross(n)/interval(j))
                    enddo
                    if(n_cross.gt.1) call ordena(cross,n_cross)
            ! Copy segments from correpondent homologues !
                    call mother(u)
                    pat=min(int(u*2+1),2)
                    sobant=0
                    do n=1,n_cross
                      genotype_t1(k,sobant+1:cross(n),j,i)=genotype(pat,sobant+1:cross(n),j,geneal(i,k))
                      sobant=cross(n)
                      select case(pat)
                        case(1)
                          pat=2
                        case(2)
                          pat=1
                      end select
                    enddo
                    genotype_t1(k,sobant+1:n_pos(j),j,i)=genotype(pat,sobant+1:n_pos(j),j,geneal(i,k))
                  enddo
                enddo
              enddo

              genotype=genotype_t1

              deallocate(genotype_t1)
              
  
    ! Cierra bucle de generaciones sin control *
      enddo   
deallocate(geneal)
        

!!!!!!Writing new genotype file and map
open(UNIT=25,FILE='base.txt')
do i=1,sum(n_ind)
  do j=1,n_chr
    write(25,'(i5,i3,1000000i2)') i,j,((genotype(k,l,j,i),k=1,2),l=1,n_pos(j))
    enddo
enddo

close(25) 

print*,'! Write to a file the new genome maps characteristics !'
open(UNIT=35,FILE='newmaps.txt')

write(35,*) '            Length (pb)  No.SNP  No.neut'
do i=1,n_chr
  write(35,'(a11,i3,i13,2i8)') 'Chromosome',i,genom_length(i),n_SNP(i),n_neut(i)
enddo
write(35,*)
do i=1,n_chr
  write(35,*) 'Chromosome',i
  do j=1,n_SNP(i)
    write(35,'(i8,i13)') SNP_pos(i,j),genom_map(i,j)
  enddo
enddo

close(35)
open(UNIT=35,FILE='newmapsneut.txt')

do i=1,n_chr
  write(35,*) 'Chromosome',i
  do j=1,n_neut(i)
    write(35,'(i8,i13)') neut_pos(i,j),genom_map_neut(i,j)
  enddo
enddo

close(35)

     


allocate(gen_1(2*sum(n_ind)), gen_2(2*sum(n_ind)))                  
allocate(frec(n_chr,maxval(n_pos)))
allocate(f_mol_neut(sum(n_ind),sum(n_ind)),f_mol_SNP(sum(n_ind),sum(n_ind)))
                

            ! Calculate different measures of diversity !
            ! in the initial population !

            call cpu_time(time1)

            ! Molecular coancestry and inbreeding !
            ! Expected and observed heterozygosity !
            ! Allelic frequencies !
            call mol_coan()

            call cpu_time(time2)
            time1=time2

            results= 0

            ! Summarises results !
            results(0,1)=sum(frec)/sum(n_pos)
            results(0,4)=(sum(frec*frec)/sum(n_pos))-(results(0,1)**2)

            ! Frequencies !
            do l=1,n_chr
                print*, 'chr frec', l
                cont=1
                do n=1,n_pos(l)
                if(n.eq.SNP_pos(l,cont)) then
	                results(0,3)=results(0,3)+frec(l,n)
	                results(0,6)=results(0,6)+frec(l,n)**2
	                cont=min(cont+1,n_SNP(l))
	            else
	                results(0,2)=results(0,2)+frec(l,n)
	                results(0,5)=results(0,5)+frec(l,n)**2
	            endif
                enddo
            enddo

            results(0,2)=results(0,2)/n_neut_total
            results(0,5)=(results(0,5)/n_neut_total)-(results(0,2)**2)


            results(0,3)=results(0,3)/n_SNP_total
            results(0,6)=(results(0,6)/n_SNP_total)-(results(0,3)**2)

            ! Expected heterozygosities !
            results(0,7)=exp_het(3,1)
            results(0,10)=exp_het(3,2)

            results(0,8)=exp_het(1,1)
            results(0,11)=exp_het(1,2)

            results(0,9)=exp_het(2,1)
            results(0,12)=exp_het(2,2)

            ! Observed heterozygosities !
            results(0,13)=obs_het(3,1)
            results(0,16)=obs_het(3,2)

            results(0,14)=obs_het(1,1)
            results(0,17)=obs_het(1,2)

            results(0,15)=obs_het(2,1)
            results(0,18)=obs_het(2,2)

            ! Segregating loci !
            results(0,19)=sum(segrega)
            results(0,20)=segrega(1)
            results(0,21)=segrega(2)

            !write(45,*)
            !write(55,*)

            ! Molecular coancestries !
            results(0,22)=sum((f_mol_neut*n_neut_total+f_mol_SNP*n_SNP_total)/(n_neut_total+n_SNP_total))
            results(0,25)=sum(((f_mol_neut*n_neut_total+f_mol_SNP*n_SNP_total)/(n_neut_total+n_SNP_total))**2)
            results(0,23)=sum(f_mol_neut)
            results(0,26)=sum(f_mol_neut*f_mol_neut)
            results(0,24)=sum(f_mol_SNP)
            results(0,27)=sum(f_mol_SNP*f_mol_SNP)

            results(0,22)=results(0,22)/(sum(n_ind)**2)
            results(0,25)=(results(0,25)/(sum(n_ind)**2))-(results(0,22)**2)
            results(0,23)=results(0,23)/(sum(n_ind)**2)
            results(0,26)=(results(0,26)/(sum(n_ind)**2))-(results(0,23)**2)
            results(0,24)=results(0,24)/(sum(n_ind)**2)
            results(0,27)=(results(0,27)/(sum(n_ind)**2))-(results(0,24)**2)
  
           ! Calculate the expected heterozygosity and the LD (at SNPs) !

cont=0
ld=0
ld2=0
dprime=0

do l=1,n_chr
  do n=1,n_SNP(l)-1
    ! Only consider segregating pairs
    if ((count(genotype(:,SNP_pos(l,n),l,:)==1).lt.2*sum(n_ind)).and. &
        (count(genotype(:,SNP_pos(l,n+1),l,:)==1).lt.2*sum(n_ind)).and. &
        (count(genotype(:,SNP_pos(l,n),l,:)==1).gt.0).and. &
        (count(genotype(:,SNP_pos(l,n+1),l,:)==1).gt.0)) then

        gen_1(1:sum(n_ind))=genotype(1,SNP_pos(l,n),l,:)
        gen_1(sum(n_ind)+1:2*sum(n_ind))=genotype(2,SNP_pos(l,n),l,:)
        gen_2(1:sum(n_ind))=genotype(1,SNP_pos(l,n+1),l,:)
        gen_2(sum(n_ind)+1:2*sum(n_ind))=genotype(2,SNP_pos(l,n+1),l,:)

        call mkdprime(dprime,r2,gen_1,gen_2,sum(n_ind)*2)

        cont=cont+1
        ld=ld+r2
        ld2=ld2+r2*r2
    endif
  enddo
enddo
if (cont.eq.0) then
    results(0,28:30)=0
else
    results(0,28)=ld/cont
    results(0,29)=(ld2/cont)-((ld/cont)**2)
    results(0,30)=cont
endif
   
print*,'! Write to a file results !'
open(UNIT=75,FILE='resultsfilrepl.txt')

    write(75,*) '     Fr_tot Fr_neu Fr_SNP Va_tot Va_neu Va_SNP'
    do i=0,0
        write(75,'(6f15.8)') (results(i,j),j=1,6)
    enddo
        write(75,*) '     He_tot Ho_tot He_neu Ho_neu He_SNP Ho_SNP'
    do i=0,0
        write(75,'(6f15.8)') (results(i,j),results(i,j+6),j=7,9)
    enddo
        write(75,*) '     Va_He_t Va_Ho_t Va_He_n Va_Ho_n Va_He_S Va_Ho_S'
    do i=0,0
        write(75,'(6f15.8)') (results(i,j),results(i,j+6),j=10,12)
    enddo
        write(75,*) '      seg_tot  seg_neu  seg_SNP'
    do i=0,0
        write(75,'(3i9)') (int(results(i,j)),j=19,21)
    enddo
        write(75,*) '      fM_tot  V_fM_t fM_neu  V_fM_n fM_SNP  V_fM_S'
    do i=0,0
        write(75,'(6f15.8)') (results(i,j),results(i,j+3),j=22,24)
    enddo
        write(75,*) '      r2_SNP  Va_r2_SNP  seg_pairs_SNP'
    do i=0,0
        write(75,'(2f9.4,i9)') results(i,28),results(i,29),int(results(i,30))
    enddo
    write(75,*) '     Fr_mafSNP  Va_mafSNP     Fr_mafneut  Va_mafneut'
        write(75,'(4f9.4)') maf,varmaf,mafneut,varmafneut
    
   close(75) 
   
deallocate(gen_1)
deallocate(gen_2)                  
deallocate(frec)
deallocate(f_mol_neut)
deallocate(f_mol_SNP)
deallocate(genotype)
deallocate(SNP_pos)
deallocate(genom_map)

diez=30
son=sum(n_ind)

!open(UNIT=55,FILE='pop_input.txt')
!    write(55,*) n_ind,son,son,n_chr,diez
!    write(55,*) n_SNP_total,n_neut_total
!    write(55,*) total_genom_length
!    write(55,*)matriz
!close(55)

end program filtrobase
    
    
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

    
    
subroutine frequency(nparent,n_chr,n_pos,n_SNP,n_neut,SNP_pos,genotype,pOBSS,pOBSSneut)
! Calculating observed  allele frequencies     
!pOBSS=observed allele frequencies

implicit none

integer :: i,j,k,l,m
integer, intent(in) ::nparent,n_chr
integer, dimension(n_chr), intent(in) :: n_SNP,n_pos,n_neut
integer, dimension(n_chr,maxval(n_SNP)), intent(in) :: SNP_pos
integer*1, dimension(2,maxval(n_pos),n_chr,nparent),  intent(in) :: genotype

integer, allocatable, dimension(:) :: contB
double precision, dimension(sum(n_SNP)),  intent(out) :: pOBSS
double precision, dimension(sum(n_neut)),  intent(out) :: pOBSSneut


allocate(contB(maxval(n_pos)))



l=0
i=0
contB=0    
do j=1,n_chr
    m=0
    do k=1,n_pos(j)
        m=m+1
        if(m.le.n_SNP(j)) then
            if(k.eq.SNP_pos(j,m)) then
                l=l+1
                contB(k)=count((genotype(:,k,j,:))==1)
                pOBSS(l)=dble(contB(k))/dble((nparent)*2)
            else
                m=m-1
                i=i+1
                contB(k)=count((genotype(:,k,j,:))==1)
                pOBSSneut(i)=dble(contB(k))/dble((nparent)*2)
            endif
        else
            i=i+1
            contB(k)=count((genotype(:,k,j,:))==1)
            pOBSSneut(i)=dble(contB(k))/dble((nparent)*2)
        endif
    end do   
end do


  

deallocate(contB)


    


  
    end  subroutine frequency
    
    
    
    

subroutine mol_coan()
use variables
implicit none
integer cont,tipo,ibs_neut,ibs_SNP,k1,k2
double precision p

frec=0

! Frequencies !
do i=1,n_chr
  do j=1,n_pos(i)
    frec(i,j)=count(genotype(:,j,i,:)==1)
  enddo
enddo
frec=frec/(2*sum(n_ind))
exp_het=0
obs_het=0

call cpu_time(time2)

time1=time2


! Expected heterozygosity for neutral, SNPs and global !
segrega=0
do i=1,n_chr
  !write(45,*) 'Chromosome',i
  !write(55,*) 'Chromosome',i
  cont=1
  do j=1,n_pos(i)
    p=frec(i,j)
    exp_het(3,1)=exp_het(3,1)+(2*p*(1-p))
    exp_het(3,2)=exp_het(3,2)+((2*p*(1-p))**2)
    obs_het(3,1)=obs_het(3,1)+(dble(count(genotype(1,j,i,:)==genotype(2,j,i,:)))/sum(n_ind))
    obs_het(3,2)=obs_het(3,2)+((dble(count(genotype(1,j,i,:)==genotype(2,j,i,:)))/sum(n_ind))**2)
	if(j.eq.SNP_pos(i,cont)) then
      exp_het(2,1)=exp_het(2,1)+(2*p*(1-p))
      exp_het(2,2)=exp_het(2,2)+((2*p*(1-p))**2)
      obs_het(2,1)=obs_het(2,1)+(dble(count(genotype(1,j,i,:)==genotype(2,j,i,:)))/sum(n_ind))
      obs_het(2,2)=obs_het(2,2)+((dble(count(genotype(1,j,i,:)==genotype(2,j,i,:)))/sum(n_ind))**2)
      cont=min(cont+1,n_SNP(i))
	  tipo=1
	  segrega(2)=segrega(2)+int(.99999+(2*p*(1-p)))
	else
      exp_het(1,1)=exp_het(1,1)+(2*p*(1-p))
      exp_het(1,2)=exp_het(1,2)+((2*p*(1-p))**2)
      obs_het(1,1)=obs_het(1,1)+(dble(count(genotype(1,j,i,:)==genotype(2,j,i,:)))/sum(n_ind))
      obs_het(1,2)=obs_het(1,2)+((dble(count(genotype(1,j,i,:)==genotype(2,j,i,:)))/sum(n_ind))**2)
	  tipo=0
	  segrega(1)=segrega(1)+int(.99999+(2*p*(1-p)))
	endif
!	write(45,'(i2,f7.4)') tipo,p
!	write(55,'(i2,2f7.4)') tipo,2*p*(1-p),1-(dble(count(genotype(1,j,i,:)==genotype(2,j,i,:)))/sum(n_ind))
  enddo
enddo



exp_het(1,1)=exp_het(1,1)/n_neut_total
exp_het(2,1)=exp_het(2,1)/n_SNP_total
exp_het(3,1)=exp_het(3,1)/sum(n_pos)

exp_het(1,2)=(exp_het(1,2)/n_neut_total)-(exp_het(1,1)**2)
exp_het(2,2)=(exp_het(2,2)/n_SNP_total)-(exp_het(2,1)**2)
exp_het(3,2)=(exp_het(3,2)/sum(n_pos))-(exp_het(3,1)**2)

obs_het(1,1)=obs_het(1,1)/n_neut_total
obs_het(2,1)=obs_het(2,1)/n_SNP_total
obs_het(3,1)=obs_het(3,1)/sum(n_pos)

obs_het(1,2)=(obs_het(1,2)/n_neut_total)-(obs_het(1,1)**2)
obs_het(2,2)=(obs_het(2,2)/n_SNP_total)-(obs_het(2,1)**2)
obs_het(3,2)=(obs_het(3,2)/sum(n_pos))-(obs_het(3,1)**2)
obs_het(:,1)=1-obs_het(:,1)

call cpu_time(time2)

time1=time2

!do
!exit
!if((i_gen.eq.n_genopt).or.(i_gen.eq.0)) then
do i=1,sum(n_ind)
  do j=i,sum(n_ind)
	ibs_neut=0
	ibs_SNP=0
    do l=1,n_chr
      cont=1
      do n=1,n_pos(l)
	    if(n.eq.SNP_pos(l,cont)) then
          do k1=1,2
            do k2=1,2
              if(genotype(k1,n,l,i).eq.genotype(k2,n,l,j)) ibs_SNP=ibs_SNP+1
            enddo
          enddo
		  cont=min(cont+1,n_SNP(l))
		else
          do k1=1,2
            do k2=1,2
              if(genotype(k1,n,l,i).eq.genotype(k2,n,l,j)) ibs_neut=ibs_neut+1
            enddo
          enddo
		endif
      enddo
    enddo
    f_mol_SNP(j,i)=dble(ibs_SNP)/(4*n_SNP_total)
    f_mol_SNP(i,j)=f_mol_SNP(j,i)
    f_mol_neut(j,i)=dble(ibs_neut)/(4*n_neut_total)
    f_mol_neut(i,j)=f_mol_neut(j,i)
  enddo
enddo
!endif
!enddo

inb_mol=0
do i=1,n_ind(1)
  inb_mol(1,1)=inb_mol(1,1)+f_mol_neut(i,i)
  inb_mol(3,1)=inb_mol(3,1)+f_mol_neut(i,i)
  inb_mol(1,2)=inb_mol(1,2)+f_mol_SNP(i,i)
  inb_mol(3,2)=inb_mol(3,2)+f_mol_SNP(i,i)
enddo
do i=n_ind(1)+1,sum(n_ind)
  inb_mol(2,1)=inb_mol(2,1)+f_mol_neut(i,i)
  inb_mol(3,1)=inb_mol(3,1)+f_mol_neut(i,i)
  inb_mol(2,2)=inb_mol(2,2)+f_mol_SNP(i,i)
  inb_mol(3,2)=inb_mol(3,2)+f_mol_SNP(i,i)
enddo
do i=1,2
  inb_mol(i,:)=inb_mol(i,:)/n_ind(i)
enddo
inb_mol(3,:)=inb_mol(3,:)/sum(n_ind)

coanc_mol(1,1)=sum(f_mol_neut(1:n_ind(1),1:n_ind(1)))/(n_ind(1)**2)
coanc_mol(2,1)=sum(f_mol_neut(n_ind(1)+1:sum(n_ind),n_ind(1)+1:sum(n_ind)))/(n_ind(2)**2)
coanc_mol(3,1)=sum(f_mol_neut)/(sum(n_ind)**2)

coanc_mol(1,2)=sum(f_mol_SNP(1:n_ind(1),1:n_ind(1)))/(n_ind(1)**2)
coanc_mol(2,2)=sum(f_mol_SNP(n_ind(1)+1:sum(n_ind),n_ind(1)+1:sum(n_ind)))/(n_ind(2)**2)
coanc_mol(3,2)=sum(f_mol_SNP)/(sum(n_ind)**2)


call cpu_time(time2)

time1=time2

    end subroutine
    
    

!-------------------
SUBROUTINE mkdprime (dprime, r2, gen_1, gen_2,n_ind)
!-------------------
! Computes D prime and R2 linkage disequilibrium measures
! between loci gen_1 and gen_2
double precision, parameter :: tol=0.0001
integer :: i, j, n_ind, gen_1, gen_2
double precision:: pi, pj, dij, dmax, dprime, r2
dimension gen_1(n_ind),gen_2(n_ind)

r2 = 0.
dprime = 0.

do i=minval(gen_1(:)), maxval(gen_1(:))
    pi = count(gen_1(:)==i)
    pi = pi/n_ind
    do j=minval(gen_2(:)), maxval(gen_2(:))
       pj = count(gen_2(:)==j)
       pj = pj / n_ind
       if (pi>tol .and. pj>tol .and. (1-pi)>tol .and. (1-pj)>tol) then
           dij = count((gen_1(:)==i) .and. (gen_2(:)==j))
           dij = dij/n_ind - pi*pj
           if (dij<0.0) then
              dmax = min (pi*pj,(1-pi)*(1-pj) )
           else
              dmax = min (pi*(1-pj), pj*(1-pi) )
           endif
           if (dmax>0.) dprime = dprime + pi*pj*abs(dij)/dmax
           r2 = r2 + dij**2 / (pi*pj)
       endif
    enddo
enddo
!--------------
    END SUBROUTINE
!--------------
     
    
    
    
      subroutine mother(u)
      Common /seeds/x(8),y(8),i16,i32
      integer x,y,carx,cary
      real u,w,t
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Randomisation of a vector !
subroutine desordena(ix,nx)
dimension ix(nx)
do i=1,nx
  call mother(u)
  n=min(int(u*nx)+1,nx)
  ic=ix(i)
  ix(i)=ix(n)
  ix(n)=ic
enddo
    end subroutine
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Sort a vector in ascending order *
subroutine ordena(ix,nx)
implicit none
integer nx,i,j,ix,ic
dimension ix(nx)
ic=0
do i=1,nx-1
  do j=i+1,nx
    if(ix(j).lt.ix(i)) then
      ic=ix(i)
      ix(i)=ix(j)
      ix(j)=ic
    endif
  enddo
enddo
    end subroutine

subroutine poisson(u,nc) !samples number of poissons from a poisson distribution 
!==================== 
INTEGER :: nc 
real :: u
double precision :: x,y 
x=EXP(-1.0d0) !variance of a poission is equal to its mean, mean is 1 recombination per Morgan 
y=x 
nc=0 
DO 
IF(x>u) exit 
nc=nc+1 
y=y/REAL(nc) 
x=x+y 
END do 
    END subroutine poisson 