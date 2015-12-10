function [cp1,cp2,varargout] = get_cp_of_planes(R1,R2,t1,t2,lim)
    xlim = lim{1}; % [ 0 W ]
    ylim = lim{2}; % [ 0 H ]

    [cp1,cp2,st] = get_cp(R1,R2,t1,t2,xlim,ylim);

    if nargout == 3
        varargout{1} = st;
    end
end

% r11 u + r12 v + t11 = r21 u' + r22 v' + t21
function [uv1,uv2] = get_uv_by_U(R1,R2,t1,t2,U)
    lhs = [ R1(:,2) -R2(:,1) -R2(:,2)];
    rhs = t2 - t1 - U*R1(:,1);
    sol = lhs\rhs;
    uv1 = [U;sol(1)];
    uv2 = sol(2:3);
end

function [uv1,uv2] = get_uv_by_V(R1,R2,t1,t2,V)
    lhs = [ R1(:,1) -R2(:,1) -R2(:,2)];
    rhs = t2 - t1 - V*R1(:,2);
    sol = lhs\rhs;
    uv1 = [sol(1);V];
    uv2 = sol(2:3);
end

function [cp1,cp2,st] = get_cp(R1,R2,t1,t2,xlim,ylim)
    st = 0;

    % specify u1 ==============================================
    xmin = xlim(1);
    xmax = xlim(2);

    ymin = ylim(1);
    ymax = ylim(2);

    % -xmin:xmax
    [uv10,uv20] = get_uv_by_U(R1,R2,t1,t2,xmin);
    [uv11,uv21] = get_uv_by_U(R1,R2,t1,t2,xmax);
    % -ymin:ymax
    [uv12,uv22] = get_uv_by_V(R1,R2,t1,t2,ymin);
    [uv13,uv23] = get_uv_by_V(R1,R2,t1,t2,ymax);

    % board1基準, 2はひとまず無視
    cp1.uv = [];
    cp2.uv = [];
    if is_valid_uv(uv10(1),uv10(2),xlim,ylim)
        cp1.uv(:,end+1) = uv10;
        cp2.uv(:,end+1) = uv20;
    end
    if is_valid_uv(uv11(1),uv11(2),xlim,ylim)
        cp1.uv(:,end+1) = uv11;
        cp2.uv(:,end+1) = uv21;
    end
    if is_valid_uv(uv12(1),uv12(2),xlim,ylim)
        cp1.uv(:,end+1) = uv12;
        cp2.uv(:,end+1) = uv22;
    end
    if is_valid_uv(uv13(1),uv13(2),xlim,ylim)
        cp1.uv(:,end+1) = uv13;
        cp2.uv(:,end+1) = uv23;
    end

    if size(cp1.uv,2) ~= 2
        st = -1
        disp '[ERROR] cannot get valid cross point'
        return
    end

    [point3ds] = Util.get_3d_points({R1,R2},{t1,t2},{cp1.uv,cp2.uv},[1 2]);
    cp1.xyz = point3ds{1};
    cp2.xyz = point3ds{2};
end

function [flag] = is_valid_uv(u, v, xlim, ylim)
    valid_u = xlim(1) <= u && u <= xlim(2);
    valid_v = ylim(1) <= v && v <= ylim(2);
    flag = valid_u && valid_v;
end
