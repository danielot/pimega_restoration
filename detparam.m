function detector = detparam(model, serial_number, assembly)

switch model
    case 'pimega_540D'
        switch serial_number
            case '2'
                switch assembly
                    otherwise
                        detector.px_array = [256 256];           % [pixels]
                        detector.chip_array = [12 12];           % [chips]
                        detector.hexa_array = [2 12];            % [hexas]
                        detector.module_array = [2 2];           % [modules]
                        
                        detector.chip_gap = 3;                   % [pixels]
                        
                        detector.module_gap_x = [                % [pixels]
                            0   8
                            4   3
                            ];
                        
                        detector.module_gap_y = [                % [pixels]
                            3   0
                            4   8
                            ];
                        
                        detector.hexa_gap = cell(detector.module_array);  % [pixels]
                        detector.hexa_gap{1,1} = [0 0 0 0 0];    % Module 1
                        detector.hexa_gap{2,1} = [0 0 0 0 0];    % Module 2
                        detector.hexa_gap{1,2} = [0 0 0 0 0];    % Module 3
                        detector.hexa_gap{2,2} = [0 0 0 0 0];    % Module 4
                        
                        detector.hexa_tilt = 6.87*pi/180;        % [radians]
                        
                        detector.pixel_size = 55e-6;             % [meters]
                end
        end
end