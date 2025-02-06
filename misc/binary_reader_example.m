% Example of using binary_reader_wrapper
% Need to change "format" to 2 when post-processing the example cases

clear; clc; close all;

mfcPath = '..';

%% 1D

binDir = fullfile(mfcPath, 'examples', '1D_acoustic_dipole', 'binary');
ti = 0;
tf = 250;
tDelta = 1;

pres = binary_reader_wrapper(binDir, ti, tf, tDelta, 1);

figure;
contourf(pres)

%% 2D

binDir  = fullfile(mfcPath, 'examples', '2D_acoustic_support5', 'binary');
ti = 0;
tf = 200;
tDelta = 10;

pres = binary_reader_wrapper(binDir, ti, tf, tDelta, 2);

figure;
for i = 1:9
    subplot(3,3,i)
    contourf(squeeze(pres(:, :, 2*i))');
end

%% 3D

binDir  = fullfile(mfcPath, 'examples', '3D_acoustic_support7', 'binary');
ti = 0;
tf = 100;
tDelta = 5;

pres = binary_reader_wrapper(binDir, ti, tf, tDelta, 3);

figure;
for i = 1:9
    subplot(3,3,i)
    contourf(squeeze(pres(:, :, 25, 2*i))');
end