    <DEM2D - simulate and compare different particle models in 2D>
    Copyright (C) <2021>  <Jonathan Jahnke>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

% Particle demonstration code with different collision solvers

% Default: DEM modell with linear spring dashpot interaction (par.algortithm = 'LIN')

% For a first start run DEM2Dmain.m
% To change the number of particles or different contact parameters, etc open DEM2Dparam.m
% To simulate the same particle setting (initial condition), set LoadData = true in DEM2Dmain.m 

% Rotations in beta stage: (some bugs pending), to turn it off set par.considerRotations = false
% Frozen particles for linear penalty based DEM
% Merged particles for linear penalty based DEM in beta stage
% Tool interaction for linear penalty based DEM

% To change the algorithm use par.algorithm = 'LIN' etc.
% Position Based Dynamics (PBD) in beta stage
% nonsmooth Projected Gauss Jacobi (PGJ)in beta stage
% nonsmooth Projected Gauss Seidel (PGS) in beta stage
% penalty based Hertz Mindlin Deresievicz DEM (HMD) in beta state