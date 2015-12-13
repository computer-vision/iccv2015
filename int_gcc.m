%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2015, Mai Nishimura, Shohei Nobuhara, and Takashi Matsuyama
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%    * Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in the
%      documentation and/or other materials provided with the distribution.
%    * Neither the name of the Graduate School of Informatics, Kyoto
%      University, Japan nor the names of its contributors may be used to
%      endorse or promote products derived from this software without specific
%      prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

% simulation data set
SIMDT = {};
SIMDT.RESOLUTION.x  = 2048;
SIMDT.RESOLUTION.y  = 1024;
SIMDT.SCALE         = 2000;

INTERSECTION.USE_POINT_NUM_ON_LINE = 10;

% generate potures for planes
[SIMDT.Rs, SIMDT.Ts] = Util.gen_postures_for_planes(SIMDT.RESOLUTION.x/SIMDT.SCALE/10);

% generate sample points
[INTERSECTION.Ps] = IntGCC.gen_samplepos_from_intersection(  SIMDT.Rs, SIMDT.Ts,...
                                                              {[-SIMDT.RESOLUTION.x/SIMDT.SCALE/2 SIMDT.RESOLUTION.x/SIMDT.SCALE/2],...
                                                              [-SIMDT.RESOLUTION.y/SIMDT.SCALE/2 +SIMDT.RESOLUTION.y/SIMDT.SCALE/2 ]}, ...
                                                              INTERSECTION.USE_POINT_NUM_ON_LINE  );

[RESULT] = IntGCC.intersection_based_calibration(INTERSECTION.Ps);

% confirm plot
subplot(1,3,1)
title('board \Phi_0')
scatter(INTERSECTION.Ps{1}.p0p1(1,:), INTERSECTION.Ps{1}.p0p1(2,:), 'r', 'filled');
hold on
scatter(INTERSECTION.Ps{1}.p2p3(1,:), INTERSECTION.Ps{1}.p2p3(2,:), 'r', 'filled');
legend('p0p1','p2p3')

subplot(1,3,2)
title('board \Phi_1')
scatter(INTERSECTION.Ps{2}.p0p1(1,:), INTERSECTION.Ps{2}.p0p1(2,:), 'g', 'filled');
hold on
scatter(INTERSECTION.Ps{2}.p4p5(1,:), INTERSECTION.Ps{2}.p4p5(2,:), 'g', 'filled');
legend('p0p1','p4p5')

subplot(1,3,3)
title('board \Phi_2')
scatter(INTERSECTION.Ps{3}.p2p3(1,:), INTERSECTION.Ps{3}.p2p3(2,:), 'b', 'filled');
hold on
scatter(INTERSECTION.Ps{3}.p4p5(1,:), INTERSECTION.Ps{3}.p4p5(2,:), 'b', 'filled');
legend('p2p3','p4p5')

figure
subplot(1,2,1)
Util.plot_rectangles_withOF(SIMDT.Rs, SIMDT.Ts, [1 2 3], SIMDT.RESOLUTION.x, SIMDT.RESOLUTION.y, 'rgb', [-SIMDT.RESOLUTION.x/2, -SIMDT.RESOLUTION.y/2])
Util.plot_rectangles_withOF({eye(3), RESULT{1}.R1.mod, RESULT{1}.R2.mod}, {zeros(3,1), RESULT{1}.t1, RESULT{1}.t2}, [1 2 3], SIMDT.RESOLUTION.x, SIMDT.RESOLUTION.y, 'kkk', [-SIMDT.RESOLUTION.x/2, -SIMDT.RESOLUTION.y/2])
subplot(1,2,2)
Util.plot_rectangles_withOF(SIMDT.Rs, SIMDT.Ts, [1 2 3], SIMDT.RESOLUTION.x, SIMDT.RESOLUTION.y, 'rgb', [-SIMDT.RESOLUTION.x/2, -SIMDT.RESOLUTION.y/2])
Util.plot_rectangles_withOF({eye(3), RESULT{2}.R1.mod, RESULT{2}.R2.mod}, {zeros(3,1), RESULT{2}.t1, RESULT{2}.t2}, [1 2 3], SIMDT.RESOLUTION.x, SIMDT.RESOLUTION.y, 'kkk', [-SIMDT.RESOLUTION.x/2, -SIMDT.RESOLUTION.y/2])
