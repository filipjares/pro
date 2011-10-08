function motion()
    %% 1.
    % Simulujte pohyb se zadanou maticí R a translací o, kde o má
    % význam o_{\beta'}, jako v rovnici 3.8 v [PRO-2011-Lecture-02.pdf].
    % Uvažujte bázi \beta jako standardní bázi.

    % matice rotace, priblizna
    R = [0.8047   -0.5059   -0.3106; ...
         0.3106    0.8047   -0.5059; ...
         0.5059    0.3106    0.8047];

    % matice rotace, lepe:
    [U,D,V] = svd(R);
    R = U*V';

    % translace;
    o = [1;1;1];
    o = -R'*o;

    %% 2.
    % Najděte souřadnice vektorů \beta' v \beta a naopak.
    % Nakreslete vektory \beta a \beta' ve standardní bázi.

    % Odpověďi:
    % - souřadnice vektorů \beta v \beta' jsou sloupce matice R
    % - souřadnice vektorů \beta' v \beta jsou sloupce transponované matice R
    %   (resp. řádky matice R)

    % pocatek a souradne soustav beta a beta' (beta2 = beta');
    beta = [ ...
        1 0 0; ...
        0 1 0; ...
        0 0 1];
    beta2 = R';
    
    figure(1); % nakreslíme souřadné soustavy ve standardní bázi:
    plot_coordinate_frame([0 0 0]', beta, 'k', 'beta', 'O'); hold on;
    plot_coordinate_frame([0 0 0]', beta2, 'b', 'beta'''); hold off;
    grid on;
    axis equal;

    %% 3.Nakreslete souřadné soustavy (O=0,\beta) a (O',\beta').
    
    figure(2);
    plot_coordinate_frame([0 0 0]', beta, 'k', 'beta', 'O'); hold on;
    plot_coordinate_frame(o, beta2, 'b', 'beta''', 'O''');
    grid on;
    axis equal;
    
    %% 4.
    % Nakreslete polohové vektory bodu X = [1;2;3]
    % vzhledem k (O=0,\beta) a (O',\beta').
    
    X = [1 2 3]';
    X2 = X - o;
    
    plot_vector([0 0 0]', X, 'g', 'x');
    plot_vector(o, X2, 'm', 'x''');
    
    %% 5.
    % Nakreslete bod Y do kterého se pohne bod X.
    
    Y = R'*X + o - o;
    
    plot_vector(o, Y, 'k', 'y''');
    hold off;
    grid on;
    axis equal;
    
end

function plot_coordinate_frame(origin, base, color, vector_label, origin_label)

    if nargin < 3
        color = 'k'
    end
    
    scale = 0;
    O = origin*[1 1 1]; % 3x origin cooridnate
    quiver3(O(1,:), O(2,:), O(3,:), base(1,:), base(2,:), base(3,:), scale, color, 'LineWidth', 2);
    
    if nargin > 3
        for i = 1:3
            text(origin(1) + 0.5*base(1,i), origin(2) + 0.5*base(2,i), origin(3) + 0.5*base(3,i),  [vector_label '_' num2str(i)]);  %,  'BackgroundColor', 'white');
        end
    end
    
    if nargin > 4
        text_pos = origin - 0.05 * (base(:,1) + base(:,2) + base(:,3));
        text(text_pos(1), text_pos(2), text_pos(3), origin_label, 'BackgroundColor', 'white');
    end
end

function plot_vector(position, vector, color, label)

    if nargin < 3
        color = 'k';
    end

    scale = 0;
    quiver3(position(1), position(2), position(3), vector(1), vector(2), vector(3), scale, color, 'LineWidth', 2);
    
    if nargin > 3
        text(position(1) + 0.5*vector(1), position(2) + 0.5*vector(2), position(3) + 0.5*vector(3), label, 'BackgroundColor', 'white');
    end

end