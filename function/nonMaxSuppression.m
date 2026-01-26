function G_nms = nonMaxSuppression(Gmag, Gdir)

[m, n] = size(Gmag);
G_nms = zeros(m,n);

% 將角度轉到 [0,180)
Gdir(Gdir < 0) = Gdir(Gdir < 0) + 180;

for i = 2:m-1
    for j = 2:n-1
        angle = Gdir(i,j);
        mag = Gmag(i,j);

        % 根據方向判斷比較像素
        if ((angle >= 0 && angle < 22.5) || (angle >= 157.5 && angle <= 180))
            q = Gmag(i, j+1);
            r = Gmag(i, j-1);
        elseif (angle >= 22.5 && angle < 67.5)
            q = Gmag(i+1, j-1);
            r = Gmag(i-1, j+1);
        elseif (angle >= 67.5 && angle < 112.5)
            q = Gmag(i+1, j);
            r = Gmag(i-1, j);
        else
            q = Gmag(i-1, j-1);
            r = Gmag(i+1, j+1);
        end

        if mag >= q && mag >= r
            G_nms(i,j) = mag;
        else
            G_nms(i,j) = 0;
        end
    end
end

end
