function R = resampleBoundary(P, N)
    d = [0; cumsum(sqrt(sum(diff(P).^2,2)))];
    d = d / d(end);
    t = linspace(0,1,N);
    Rx = interp1(d, P(:,1), t, 'linear');
    Ry = interp1(d, P(:,2), t, 'linear');
    R = [Rx' Ry'];
end