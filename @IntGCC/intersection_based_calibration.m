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

function [result] = intersection_based_calibration(POSs, varargin)
  % POSs => { struct('p0p1', [2x10], 'p2p3', [2x10]), ... }
  ptnum = size(POSs{1}.p0p1,2);

  msub1 = [];
  msub2 = [];
  for idx=1:2:ptnum-1
    ptidx = [idx idx+1];
    [m1,m2] = get_mat_from_pos(POSs, ptidx);

    m1row = size(m1,1);
    m2row = size(m2,1);
    msub1(end+1:end+m1row ,:) = m1;
    msub2(end+1:end+m2row ,:) = m2;
  end

  % msub-inを追加するよ
  msub1(end+1:end+2, :) = get_inner_product_row(POSs);

  % Mx = b
  Msub1 = msub1(:,1:end-1); % 14x12 matrix RANK = 12
  b = -msub1(:,end);
  % Supecial Solution
%     Xsub1 = inv(Msub1'*Msub1)*Msub1'*b;
  Xsub1 = Msub1\b;

  % Mx = 0
  Msub2 = msub2(:,1:end-1); % RANK = 5
  [U,S,V] = svd(Msub2);
  Xsub2 = V(:,end);

%     Xsub2 = Msub2\b;
  result = gen_RT_from_2x(Xsub1,Xsub2);
end


% PCA
function [ dir_vec_origin, dir_vec ] = get_2dvector_from_points(points)
  assert(size(points,1) > size(points,2));
  % points = [ u v ; u v ; u v]
  covmat = cov(points);
  [V,D] = eig(covmat);

  dir_vec_origin = mean(points)';
  dir_vec = V(:,2); % first eigen vector
  dir_vec = dir_vec/norm(dir_vec); % constrains to be unit vector
end

function [ msub1, msub2 ] = get_mat_from_pos(POSs,idx)
  pos0 = POSs{1};
  pos1 = POSs{2};
  pos2 = POSs{3};

  % board 0
  u00 = pos0.p0p1(1,idx(1));
  v00 = pos0.p0p1(2,idx(1));
  u01 = pos0.p0p1(1,idx(2));
  v01 = pos0.p0p1(2,idx(2));

  u02 = pos0.p2p3(1,idx(1));
  v02 = pos0.p2p3(2,idx(1));
  u03 = pos0.p2p3(1,idx(2));
  v03 = pos0.p2p3(2,idx(2));

  % board 1
  u10 = pos1.p0p1(1,idx(1));
  v10 = pos1.p0p1(2,idx(1));
  u11 = pos1.p0p1(1,idx(2));
  v11 = pos1.p0p1(2,idx(2));

  u12 = pos1.p4p5(1,idx(1));
  v12 = pos1.p4p5(2,idx(1));
  u13 = pos1.p4p5(1,idx(2));
  v13 = pos1.p4p5(2,idx(2));

  % board 2
  u20 = pos2.p2p3(1,idx(1));
  v20 = pos2.p2p3(2,idx(1));
  u21 = pos2.p2p3(1,idx(2));
  v21 = pos2.p2p3(2,idx(2));

  u22 = pos2.p4p5(1,idx(1));
  v22 = pos2.p4p5(2,idx(1));
  u23 = pos2.p4p5(1,idx(2));
  v23 = pos2.p4p5(2,idx(2));

  % Mx = b
  msub1 = [
              u10,v10,0,0,0,0,0,0,1,0,0,0,-u00;
              0,0,u10,v10,0,0,0,0,0,1,0,0,-v00;
              u11,v11,0,0,0,0,0,0,1,0,0,0,-u01;
              0,0,u11,v11,0,0,0,0,0,1,0,0,-v01;
              0,0,0,0,u20,v20,0,0,0,0,1,0,-u02;
              0,0,0,0,0,0,u20,v20,0,0,0,1,-v02;
              0,0,0,0,u21,v21,0,0,0,0,1,0,-u03;
              0,0,0,0,0,0,u21,v21,0,0,0,1,-v03;
              u12,v12,0,0,-u22,-v22,0,0,1,0,-1,0,0;
              0,0,u12,v12,0,0,-u22,-v22,0,1,0,-1,0;
              u13,v13,0,0,-u23,-v23,0,0,1,0,-1,0,0;
              0,0,u13,v13,0,0,-u23,-v23,0,1,0,-1,0
              ];

  % Mx = 0
  msub2 = [   u10,v10,0,0,1,0,0;
              u11,v11,0,0,1,0,0;
              0,0,u20,v20,0,1,0;
              0,0,u21,v21,0,1,0;
              u12,v12,-u22,-v22,1,-1,0;
              u13,v13,-u23,-v23,1,-1,0 ];
end

function [msub1_in] = get_inner_product_row(POSs)
  pos0 = POSs{1};
  pos1 = POSs{2};
  pos2 = POSs{3};

  % board 0
  uv01 = pos0.p0p1;
  uv02 = pos0.p2p3;

  % board 1
  uv10 = pos1.p0p1;
  uv12 = pos1.p4p5;

  % board 2
  uv20 = pos2.p2p3;
  uv21 = pos2.p4p5;

  % get unit vector from 2 points
  [ center, UV01 ] = get_2dvector_from_points(uv01');
  [ center, UV02 ] = get_2dvector_from_points(uv02');
  [ center, UV10 ] = get_2dvector_from_points(uv10');
  [ center, UV12 ] = get_2dvector_from_points(uv12');
  [ center, UV20 ] = get_2dvector_from_points(uv20');
  [ center, UV21 ] = get_2dvector_from_points(uv21');

  % board 0
  u00 = 0;
  v00 = 0;
  u01 = UV01(1);
  v01 = UV01(2);

  u02 = 0;
  v02 = 0;
  u03 = UV02(1);
  v03 = UV02(2);

  % board 1
  u10 = 0;
  v10 = 0;
  u11 = UV10(1);
  v11 = UV10(2);

  u12 = 0;
  v12 = 0;
  u13 = UV12(1);
  v13 = UV12(2);

  % board 2
  u20 = 0;
  v20 = 0;
  u21 = UV20(1);
  v21 = UV20(2);

  u22 = 0;
  v22 = 0;
  u23 = UV21(1);
  v23 = UV21(2);

  msub1_in = [
              u01*u13-u00*u13-u01*u12+u00*u12,u01*v13-u00*v13-u01*v12+u00*v12,u13*v01-u12*v01-u13*v00+u12*v00,v01*v13-v00*v13-v01*v12+v00*v12,0,0,0,0,0,0,0,0,-v11*v13+v10*v13+v11*v12-v10*v12-u11*u13+u10*u13+u11*u12-u10*u12;
              0,0,0,0,u03*u23-u02*u23-u03*u22+u02*u22,u03*v23-u02*v23-u03*v22+u02*v22,u23*v03-u22*v03-u23*v02+u22*v02,v03*v23-v02*v23-v03*v22+v02*v22,0,0,0,0,-v21*v23+v20*v23+v21*v22-v20*v22-u21*u23+u20*u23+u21*u22-u20*u22
              ];
end

function [ESTIM_DATA] = gen_RT_from_2x(x_xy,x_z)
  % x_xy : [r11,r12,r14,r15,r21,r22,r24,r25,t11,t12,t21,t22]
  % x_z  : [r17,r18,r27,r28,t13,t23]
  r11 = x_xy(1);
  r12 = x_xy(2);
  r14 = x_xy(3);
  r15 = x_xy(4);
  r21 = x_xy(5);
  r22 = x_xy(6);
  r24 = x_xy(7);
  r25 = x_xy(8);
  t11 = x_xy(9);
  t12 = x_xy(10);
  t21 = x_xy(11);
  t22 = x_xy(12);

  r17 = x_z(1);
  r18 = x_z(2);
  r27 = x_z(3);
  r28 = x_z(4);
  t13 = x_z(5);
  t23 = x_z(6);

  % |r11| = 1
  x0 = real( sqrt( ( 1 - (r11^2 + r14^2) ) / r17^2 ) );
  xdata.x_xy = x_xy;
  xdata.x_z  = x_z;
  f = orthogonality_constraints(x0,xdata);
  ydata = zeros(length(f),1);
  opt = optimset('algorithm', 'levenberg-marquardt', 'TolFun', 1e-20, 'MaxFunEvals', 500000, 'MaxIter', 100000, 'TolX', 1e-20);
      [opt_x, resnorm, residual, exitflag, output] = lsqcurvefit(@orthogonality_constraints,x0,xdata,ydata,[],[],opt);

  x_z = opt_x*x_z;
  r17 = x_z(1);
  r18 = x_z(2);
  r27 = x_z(3);
  r28 = x_z(4);
  t13 = x_z(5);
  t23 = x_z(6);

  ESTIM_DATA = {};
  % x_xy + Bx_z
  ESTIM_DATA{1}.R1 = Util.gen_R_from_partial(r11,r12,r14,r15,r17,r18);
  ESTIM_DATA{1}.R2 = Util.gen_R_from_partial(r21,r22,r24,r25,r27,r28);
  ESTIM_DATA{1}.t1 = [t11;t12;t13];
  ESTIM_DATA{1}.t2 = [t21;t22;t23];
  % x_xy - Bx_z
  ESTIM_DATA{2}.R1 = Util.gen_R_from_partial(r11,r12,r14,r15,-r17,-r18);
  ESTIM_DATA{2}.R2 = Util.gen_R_from_partial(r21,r22,r24,r25,-r27,-r28);
  ESTIM_DATA{2}.t1 = [t11;t12;-t13];
  ESTIM_DATA{2}.t2 = [t21;t22;-t23];

  % Ax = b 側は反転ないはず
%     % -x_xy + Bx_z
%     ESTIM_DATA{3}.R1 = Util.gen_R_from_partial(-r11,-r12,-r14,-r15,B*r17,B*r18);
%     ESTIM_DATA{3}.R2 = Util.gen_R_from_partial(-r21,-r22,-r24,-r25,B*r27,B*r28);
%     ESTIM_DATA{3}.t1 = [-t11;-t12;B*t13];
%     ESTIM_DATA{3}.t2 = [-t21;-t22;B*t23];
%
%     % -x_xy - Bx_z
%     ESTIM_DATA{4}.R1 = Util.gen_R_from_partial(-r11,-r12,-r14,-r15,-B*r17,-B*r18);
%     ESTIM_DATA{4}.R2 = Util.gen_R_from_partial(-r21,-r22,-r24,-r25,-B*r27,-B*r28);
%     ESTIM_DATA{4}.t1 = [-t11;-t12;-B*t13];
%     ESTIM_DATA{4}.t2 = [-t21;-t22;-B*t23];
end

function f = orthogonality_constraints(x, xdata)
  x_xy = xdata.x_xy;
  x_z  = xdata.x_z;

  r11 = x_xy(1);
  r12 = x_xy(2);
  r14 = x_xy(3);
  r15 = x_xy(4);
  r21 = x_xy(5);
  r22 = x_xy(6);
  r24 = x_xy(7);
  r25 = x_xy(8);
  t11 = x_xy(9);
  t12 = x_xy(10);
  t21 = x_xy(11);
  t22 = x_xy(12);

  r17 = x_z(1);
  r18 = x_z(2);
  r27 = x_z(3);
  r28 = x_z(4);
  t13 = x_z(5);
  t23 = x_z(6);

  % |r11| = 1
  EQ1 = r11^2 + r14^2 + (x*r17)^2 - 1;
  % |r11| = 1
  EQ2 = r12^2 + r15^2 + (x*r18)^2 - 1;
  % r11.r12 = 0
  EQ3 = dot([r11 r14 x*r17], [r12 r15 x*r18]);

  % |r21| = 1
  EQ4 = r21^2 + r24^2 + (x*r27)^2 - 1;
  % |r21| = 1
  EQ5 = r22^2 + r25^2 + (x*r28)^2 - 1;
  % r21.r22 = 0
  EQ6 = dot([r21 r24 x*r27], [r22 r25 x*r28]);

  f = [EQ1;EQ2;EQ3;EQ4;EQ5;EQ6];
end
