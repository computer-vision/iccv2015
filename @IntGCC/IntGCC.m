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

% Generalized Camera Calibartion from Intersections
classdef IntGCC < handle
  properties

  end

  methods
    function obj = GCM(varargin) % CONSTRUCTOR
    end
  end

  methods(Static)

    function [POSs] = gen_samplepos_from_intersection(Rs,Ts,lim,samplenum)
      R0 = Rs{1};
      R1 = Rs{2};
      R2 = Rs{3};
      t0 = Ts{1};
      t1 = Ts{2};
      t2 = Ts{3};

      [cp11 cp21, st1] = Util.get_cp_of_planes(R1,R2,t1,t2,lim);
      [cp00 cp10, st2] = Util.get_cp_of_planes(R0,R1,t0,t1,lim);
      [cp01 cp20, st3] = Util.get_cp_of_planes(R0,R2,t0,t2,lim);

      POSs = {};

      if st1 ~= 0 || st2 ~= 0 || st3 ~= 0
        struct('st1', st1, 'st2', st2, 'st3', st3)
        return
      end

      % BOARD 0
      POSs{1}.p0p1 = Util.interpolate_2point(cp00.uv, samplenum);
      POSs{1}.p2p3 = Util.interpolate_2point(cp01.uv, samplenum);

      % BOARD 1
      POSs{2}.p0p1 = Util.interpolate_2point(cp10.uv, samplenum);
      POSs{2}.p4p5 = Util.interpolate_2point(cp11.uv, samplenum);

      % BOARD 2
      POSs{3}.p2p3 = Util.interpolate_2point(cp20.uv, samplenum);
      POSs{3}.p4p5 = Util.interpolate_2point(cp21.uv, samplenum);
    end

    % FIXME
    [result] = intersection_based_calibration(POSs, varargin)
  end
end
