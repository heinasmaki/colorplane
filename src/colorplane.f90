!
!   Colorplane: simulated annealing of a modified Potts model
!
!   Copyright (C) 2018 Sami Heinäsmäki
!
!   Colorplane is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Colorplane is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with Colorplane.  If not, see <https://www.gnu.org/licenses/>.
!
program colorplane
  use colorplane_consts
  use colorplane_utils
  use colorplane_reliability
  use colorplane_parallel
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  type(basedata) :: lattice
  type(parallel_plan) :: plan
  integer(ip), allocatable, target :: spins(:,:)
  integer(ip), pointer :: flat(:)
  integer :: device
  logical :: parallel_mode
  integer, allocatable :: sites(:)
  integer :: i,j,k,ct,cx,cy,id,is,js,pos,site,metrics,history,interior_sites
  integer(int64) :: attempts,accepted,changed,neutral
  real(dp) :: T,e_old,e_new,r,p,err,mean,interior,zero_fraction,maximum,denom
  integer(ip) :: s_old,s_new,s_ref,s_out
  real(dp) :: start,finish
  character(len=1024) :: arg
  logical :: state_exists

  call get_command_argument(1,arg)
  if (arg=='--help') then
     call usage()
     stop
  end if
  call collect_input(lattice)
  call read_options(lattice)
  call validate_input(lattice)
  call configure_backend(lattice,device)
  parallel_mode=lattice%backend/='serial'
  call check_if_restart(lattice)
  inquire(file=trim(lattice%outfil)//'.state',exist=state_exists)
  if (state_exists .and. .not. lattice%restart) error stop 'Incomplete checkpoint: .lat/.sum missing'
  if (lattice%restart) call read_checkpoint(lattice)
  call validate_input(lattice)
  call finalize_input(lattice)
  allocate(spins(lattice%Ny+2*lattice%d,lattice%Nx+2*lattice%d))
  if (.not. parallel_mode) allocate(sites(lattice%Nx*lattice%Ny))
  spins=0
  if (lattice%restart) then
     call load_lattice(trim(lattice%outfil)//'.lat',spins,lattice)
     if (lattice_hash(spins,lattice)/=lattice%checkpoint_hash) error stop 'Checkpoint lattice checksum mismatch'
     lattice%outfil=trim(lattice%outfil)//'-res'
     call check_output_root(lattice)
  else
     call check_output_root(lattice)
     call initialize_rng(lattice)
     if (len_trim(lattice%initial_file)>0) then
        call load_lattice(trim(lattice%initial_file),spins,lattice)
     else
        call create_initial_spins(lattice,spins)
     end if
  end if
  flat(1:size(spins))=>spins
  if (parallel_mode) then
     call make_plan(lattice,plan)
     lattice%tile_size=plan%tile
     lattice%color_groups=plan%groups
     call parallel_begin(flat,plan)
     print '(a,i0,a,i0)', 'Independent tiles: side=',plan%tile,', color groups=',plan%groups
  end if
  ! Record the generator state immediately before the first Monte Carlo move.
  call save_rng(lattice)
  call save_lattice(trim(lattice%outfil)//'.initial.lat',spins,lattice)
  call write_metadata(lattice)
  if (.not. parallel_mode) sites=[(i,i=1,size(sites))]
  open(newunit=metrics,file=trim(lattice%outfil)//'.metrics.csv',status='replace')
  write(metrics,'(a)') 'step,temperature,energy,interior_energy,interior_sites,zero_fraction,max_local_energy,'// &
       'acceptance_rate,energy_change_acceptance_rate,zero_energy_acceptance_rate'
  open(newunit=history,file=trim(lattice%outfil)//'.dat',status='replace')
  if (parallel_mode) then
     write(history,'(a)') '# Parallel mode: fixed-configuration samples are in .metrics.csv; no per-sweep transfers'
  else
     write(history,'(a)') '# Sweep-time pre-update energy estimate; use .metrics.csv for fixed configurations'
  end if
  start=wall_seconds()
  do k=1,lattice%tmax
     T=lattice%tvec(k)
     attempts=0
     accepted=0
     changed=0
     neutral=0
     if (parallel_mode) call parallel_reset(plan)
     do ct=1,lattice%ctmax-lattice%offset*(lattice%completed_steps+k)
        if (parallel_mode) then
           if (lattice%random_order) call shuffle_sites(plan%order)
           call parallel_sweep(flat,plan,lattice,T)
           lattice%completed_sweeps=lattice%completed_sweeps+1
           cycle
        end if
        if (lattice%random_order) call shuffle_sites(sites)
        err=0
        do pos=1,size(sites)
           site=sites(pos)
           i=1+(site-1)/lattice%Ny
           j=1+mod(site-1,lattice%Ny)
           cx=i+lattice%d
           cy=j+lattice%d
           s_old=spins(cy,cx)
           s_new=newspin(s_old,lattice%q)
           e_old=0
           e_new=0
           do id=1,lattice%dnum
              is=i+lattice%xind(id)-1
              js=j+lattice%yind(id)-1
              s_ref=spins(js,is)
              if (s_ref==s_old) e_old=e_old+lattice%dval(id)
              if (s_ref==s_new) e_new=e_new+lattice%dval(id)
           end do
           s_out=s_old
           if (e_new<e_old) then
              s_out=s_new
           else
              p=exp(-(e_new-e_old)/T)
              call random_number(r)
              if (r<p) s_out=s_new
           end if
           attempts=attempts+1
           if (s_out/=s_old) then
              accepted=accepted+1
              if (e_new==e_old) then
                 neutral=neutral+1
              else
                 changed=changed+1
              end if
           end if
           spins(cy,cx)=s_out
           err=err+e_old
        end do
        lattice%completed_sweeps=lattice%completed_sweeps+1
        write(history,'(es25.17)') err/(real(size(sites),dp)*sum(lattice%dval))
     end do
     if (parallel_mode) then
        call parallel_measure(flat,plan,lattice,mean,interior,zero_fraction,maximum,interior_sites,accepted,changed)
        attempts=int(lattice%Nx,int64)*lattice%Ny* &
             (lattice%ctmax-lattice%offset*(lattice%completed_steps+k))
        neutral=accepted-changed
     else
        call measure(spins,lattice,mean,interior,zero_fraction,maximum,interior_sites)
     end if
     lattice%errvec(k)=mean
     denom=real(attempts,dp)
     write(metrics,'(i0,3(a,es25.17),a,i0,5(a,es25.17))') &
          lattice%completed_steps+k,',',T,',',mean,',',interior, &
          ',',interior_sites,',',zero_fraction,',',maximum,',',accepted/denom,',',changed/denom,',',neutral/denom
     write(*,'(a,i0,a,es16.8)') 'Stage ',lattice%completed_steps+k,': exact energy = ',mean
  end do
  close(metrics)
  close(history)
  if (parallel_mode) call parallel_finish(flat,plan)
  finish=wall_seconds()
  call write_summary(lattice,lattice%tmax,finish-start)
  call save_lattice(trim(lattice%outfil)//'.lat',spins,lattice)
  if (lattice%writeneigh) then
     if (parallel_mode) then
        call parallel_map(plan,lattice)
     else
        call print_neighs(spins,lattice)
     end if
  end if
  call write_checkpoint(lattice,spins)
  print '(a,f16.6)', 'Annealing wall seconds: ',finish-start
end program
