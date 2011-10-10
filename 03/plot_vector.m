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