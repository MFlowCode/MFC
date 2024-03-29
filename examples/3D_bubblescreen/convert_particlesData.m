% Post-processing the particle file
% Author: Kazuki Maeda
% 09/01/2018

close all
clear all
clc

crd=pwd();

mkdir('m_data_particles');

% Time stepping information follows .py input files 
%%%%%%%%%%%%USER INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%
t_step_start=0;
t_step_save=50;
t_step_stop=3000;

%%%%%%%%%%%%USER INPUTS%%%%%%%%%%%%%%%%%%%%%%%%

j=1;
for i=t_step_start:t_step_save:t_step_stop
      filename = sprintf('./particles_data/particles_data%d.dat',i)
      tmp=sortrows(load(filename));
	  data(j).part(:,1:4)=tmp(:,1:4);
	  data(j).part(:,5:6)=tmp(:,11:12);
	  data(j).part(:,7)=tmp(:,22);
      j=j+1;
end

% data(i).part(j,k) includes data of j-th particle in i-th output file.
% k=1 : ID
% k=2:4 : spatial coordinate
% k=5 : radius
% k=6 : radial velocity
% k=7 : time

cd('./m_data_particles');

filename=strcat('particles_data.mat');
save('-V7.3',filename,'data');
cd('../');

clear data;

cd(crd);

