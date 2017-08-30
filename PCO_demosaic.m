function Ideconv = PCO_demosaic(img_rgb,coupling_matrix)

    Imeasurement = zeros(1,size(img_rgb,1)*size(img_rgb,2)/4);
    for bayerpatternIdx = 1:4
            switch bayerpatternIdx
                case 1
                    img_sub = img_rgb(1:2:end,1:2:end);
                case 2
                    img_sub = img_rgb(1:2:end,2:2:end);
                case 3
                    img_sub = img_rgb(2:2:end,1:2:end);
                case 4
                    img_sub = img_rgb(2:2:end,2:2:end);                
            end
            Imeasurement(bayerpatternIdx,:) = img_sub(:)';
    end
    Ideconv = mean(coupling_matrix,3)\Imeasurement;
    Ideconv = reshape(permute(Ideconv,[2 1]),[size(img_rgb,1)/2 size(img_rgb,2)/2 3]);
    % Rescale back to max amplitude
    Ideconv = Ideconv/max(Ideconv(:)) * max(img_rgb(:));
end