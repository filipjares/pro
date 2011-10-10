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