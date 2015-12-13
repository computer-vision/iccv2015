classdef Util < handle % need handle? probably, no.
  properties
  end
  methods
    function obj = Util(varargin) % CONSTRUCTOR
    end
  end

  methods(Static)
    function [Rs,Ts] = gen_postures_for_planes(t_range)
      Rs = {};
      Rs{1} = eye(3);
      Ts{1} = zeros(3,1);

      % T_LIMIT = SIMDT.RESOLUTION/20
      Rs{2} = Util.gen_R([Util.nprand(pi/3), Util.nprand(pi/3), Util.nprand(pi/3)]);
      Ts{2} = [Util.nprand(t_range), Util.nprand(t_range), Util.nprand(t_range)]';
      Rs{3} = Util.gen_R([Util.nprand(pi/3), Util.nprand(pi/3),Util.nprand(pi/3)]);
      Ts{3} = [Util.nprand(t_range), Util.nprand(t_range), Util.nprand(t_range)]';
    end

    function [R] = gen_R(theta)
      Rx = [  1,                0,                  0;
      0,                cos(theta(1,1)),    sin(theta(1,1));
      0,                -sin(theta(1,1)),   cos(theta(1,1))];
      Ry = [  cos(theta(1,2)),  0,                  -sin(theta(1,2));
      0,                1,                  0;
      sin(theta(1,2)),  0,                  cos(theta(1,2))];
      Rz = [  cos(theta(1,3)),  sin(theta(1,3)),    0;
      -sin(theta(1,3)), cos(theta(1,3)),    0;
      0,                0,                  1 ];
      R = Rx*Ry*Rz;
    end

    function [] = plot_rectangles_withOF(Rs, Ts, set_idx_list, w, h, colorstr, offset)
      offset = offset(:);
      base = [0   0   w   w;
      0   h   h   0];

      base(1,:) = base(1,:) + offset(1);
      base(2,:) = base(2,:) + offset(2);

      cstr_idx = 1;
      for set_idx=set_idx_list
                R = Rs{set_idx};
                T = Ts{set_idx};
                H = [R(:,1:2) T];

                base(3,:) = 1;

                points3d = H*base;
                plot3([points3d(1,:) points3d(1,1)], [points3d(2,:) points3d(2,1)],  [points3d(3,:) points3d(3,1)], colorstr(cstr_idx), 'LineWidth', 2);
                hold on
                cstr_idx = cstr_idx + 1;
      end
    end

    % FIXME
    function [out] = nprand(varargin)
      if nargin == 1
        num = varargin{1};
        out = rand*num - num/2;
      elseif nargin == 3
        num = varargin{1};
        row = varargin{2};
        col = varargin{3};
        out = num*rand(row,col) - num/2*ones(row,col);
      else
        assert(1, 'Usage: nprand(num) or nprand(num,row,col)');
      end
    end

    % FIXME
    [result] = intersection_based_calibration(POSs, varargin)

    % FIXME
    [cp1,cp2,varargout] = get_cp_of_planes(R1,R2,t1,t2,lim);

    % FIXME
    % 3d points cell
    function [points3ds] = get_3d_points(Rs,Ts,POSs,set_idx_list)
      points3ds = {}; %  共通ディスプレイポイント
      for set_idx=set_idx_list
        R = Rs{set_idx};
        T = Ts{set_idx};

        points2d = POSs{set_idx};
        points2d(3,:) = 0;

        points3d = R*points2d + repmat(T,1,size(points2d,2));

        points3ds{set_idx} = points3d;
        % points3ds = [ points3ds { points3d } ];
      end
    end

    % FIXME
    % [u1 u2; v1 v2]
    function [points] = interpolate_2point(uv, num)
      uv1 = uv(:,1);
      uv2 = uv(:,2);

      points = [];
      for iter=1:num
        points(:,end+1) = (iter*uv1 + (num-iter)*uv2)/num;
      end
    end

    function result = gen_R_from_partial(r1,r2,r4,r5,r7,r8)
      R = [ r1 r2;
      r4 r5;
      r7 r8 ];
      R = [R, cross(R(:,1),R(:,2))];
      % correct
      [U,S,V] = svd(R);
      Rmod = U*V';

      result.org = R;
      result.mod = Rmod;
    end
  end
end
