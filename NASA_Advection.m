function NASA_Advection
% Author: James Jalal Behzadi
% UNSW (university of New South Wales), Sydney, Australia
% Date: 13.12.2012
% Version: 1.0.0
clc

CFL = .4;
Nexact = 2000;
N = 100;

tEnd = 10;
xLow = -1;
xHigh = 1;
delX = (xHigh - xLow) / N;
delX2 = (xHigh - xLow) / Nexact;
[x c xExact cExact] = IVP(xLow, xHigh, N, Nexact, delX, delX2);

u = -1;
Nmove = CFL * Nexact / N;
delT = CFL * delX / abs(u);
NdelT = (tEnd-0) / delT;

% initialising video parameters
videoName = 'videosNASA/';
mkdir(videoName);
videoName = [videoName 'NASAadvection.avi'];
vidObj = VideoWriter(videoName);
vidObj.FrameRate = 10;
open(vidObj);

t = 0;
figurePlotter(t, x, c, xExact, cExact, N, CFL, Nmove);
currFrame = getframe(gcf);
writeVideo(vidObj,currFrame);
for i = 1:NdelT
    t = t + delT;
%     c = solve1stOrder(c, delT, N, delX, u);
    c = TVD_RK3(c, N, delX, u, delT);
    cExact = moveNASAwave(Nmove, xExact, cExact, Nexact, u);
    figurePlotter(t, x, c, xExact, cExact, N, CFL, Nmove);
    pause(1e-3)
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
end
close(vidObj);
end

function figurePlotter(t, x, c, xExact, cExact, N, CFL, Nmove)
plot(xExact, cExact, 'r', x, c, 'bo', ...
    'LineWidth', 2, ...
    'MarkerSize', 5)
xlim([-1 1])
ylim([-.3 1.3])
tChar = sprintf('%.4g', t);
cflChar = sprintf('%.4g', CFL);
NChar = sprintf('%.4g', N);
NmoveChar = sprintf('%.4g', Nmove);
title (['CFL = ', cflChar, '        Nmove = ', NmoveChar, '        N = ', NChar, '       t = ' , tChar, ' s'], ...
        'FontWeight', 'bold', ...
    'FontSize', 15)
end

function [xOut cOut xExact cExact] = IVP(xLow, xHigh, N, Nexact, delX, delX2)

xOut = linspace(xLow - 3*delX/2, xHigh + 3*delX/2, N+4);
cOut = qinit(xOut);

xExact = linspace(xLow + delX2/2, xHigh - delX2/2, Nexact);
cExact = qinit(xExact);
end

function cOut = qinit(xIn)

N = length(xIn);
a = .5;
Z = -.7;
delta = .005;
alpha = 10;
beta = log10(2) / (36 * delta ^ 2);
cOut = NaN * ones(N, 1);
for i = 1:N
    x = xIn(i);
    if (x >= -.8 & x <= -.6)
        cOut(i) = (G(x, beta, Z-delta) + G(x, beta, Z+delta) + 4 * G(x, beta, Z)) / 6;
    elseif (x >= -.4 & x <= -.2)
        cOut(i) = 1;
    elseif (x >= 0 & x <= .2)
        cOut(i) = 1 - abs(10 * (x - .1));
    elseif (x >= .4 & x <= .6)
        cOut(i) = (F(x, alpha, a-delta) + F(x, alpha, a+delta) + 4 * F(x, alpha, a)) / 6;
    else
        cOut(i) = 0;
    end
end

% figure(1)
% plot(xIn, cOut, 'r')
% ylim([-.3 1.3])
% close(1)
end

function cExact = moveNASAwave(Nmove, xExact, cExact, Nexact, u)

cTmp = NaN * ones(Nexact, 1);

if (u > 0)
    cTmp(1:Nmove) = cExact(Nexact - Nmove + 1:Nexact);
    cTmp(Nmove+1:Nexact) = cExact(1:Nexact - Nmove);
    cExact = cTmp;
else
    cTmp(1:Nexact-Nmove) = cExact(Nmove+1:Nexact);
    cTmp(Nexact-Nmove+1:Nexact) = cExact(1:Nmove);
    cExact = cTmp;
end

% figure(1)
% plot(xExact, cExact, 'r')
% ylim([-.3 1.3])
% close(1)
end

function out = G(x, beta, Z)
out = exp(-beta * (x-Z) ^ 2);
end

function out = F(x, alpha, a)
out = (max((1-alpha ^ 2 *(x-a) ^ 2), 0)) ^ .5;
end

function c = bc(c, N)

% left boundary
c(1) = c(N+1);
c(2) = c(N+2);
% right boundary
c(N+4) = c(4);
c(N+3) = c(3);
end

function out = psi(r)
% vanLeer's limiter
out = (r + abs(r)) / nonZeroDenom(1 + r);
end

function out = L(c, N, delX, u)
% TVD operator for spatial discretisation
out = NaN * ones(N, 1);
for i = 3:N+2
    rp1 = (c(i+1) - c(i)) / nonZeroDenom(c(i) - c(i-1));
%     rm1 = (c(i-1) - c(i-2)) / nonZeroDenom(c(i) - c(i-1));
    rm1 = (c(i) - c(i-1)) / nonZeroDenom(c(i+1) - c(i));
    rp2 = (c(i) - c(i-1)) / nonZeroDenom(c(i-1) - c(i-2));
    rm2 = (c(i+1) - c(i)) / nonZeroDenom(c(i+2) - c(i+1));
    uplus = max(u, 0);
    uminus = min(u, 0);
    out(i-2) = ...
        -1 / delX * (1 + .5 * psi(rp1) - .5 * psi(rp2) / nonZeroDenom(rp2)) * uplus * (c(i) - c(i-1)) ...
        -1 / delX * (1 + .5 * psi(rm1) - .5 * psi(rm2) / nonZeroDenom(rm2)) * uminus * (c(i+1) - c(i));
end
end

function out = nonZeroDenom(denom)
if (abs(denom) < eps)
    if (denom < 0)
        out = -eps;
    else
        out = eps;
    end
else
    out = denom;
end
end

function c = TVD_RK3(c, N, delX, u, delT)
% 3rd order TVD Runge-Kutta method

% applying boundary conditions before advancing the solution
c = bc(c, N);
cTmp1 = NaN * c;
cTmp1(3:N+2) = c(3:N+2) + delT * L(c, N, delX, u);

cTmp1 = bc(cTmp1, N);
cTmp2 = NaN * c;
cTmp2(3:N+2) = 3 / 4 * c(3:N+2) + 1 / 4 * cTmp1(3:N+2) + 1 / 4 * delT * L(cTmp1, N, delX, u);

cTmp2 = bc(cTmp2, N);
c(3:N+2) = 1 / 3 * c(3:N+2) + 2 / 3 * cTmp2(3:N+2) + 2 / 3 * delT * L(cTmp2, N, delX, u);
end

function c = solve1stOrder(c, delT, N, delX, u)
% 1st order time integration

% applying boundary conditions before advancing the solution
c = bc(c, N);
RHS= L(c, N, delX, u);
for i = 1:N
    c(i+2) = c(i+2) + delT * RHS(i);
end
end