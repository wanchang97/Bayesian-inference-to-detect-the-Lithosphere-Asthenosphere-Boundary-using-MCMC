function obs = getSolutionFromTable(fwd,m)
% m dim: np x 1
np = size(m,1);
switch np
    case 2
        found = false;
        for I1 = 1:length(fwd.table.grid)-1
            if fwd.table.grid(I1) <= m(1) && m(1) <= fwd.table.grid(I1+1)
                found = true;
                break
            end
        end
        assert(found,'m1 not found')

        found = false;
        for I2 = 1:length(fwd.table.grid)-1
            if fwd.table.grid(I2) <= m(2) && m(2) <= fwd.table.grid(I2+1)
                found = true;
                break
            end
        end
        assert(found,'m2 not found')

        obs = fwd.table.solVector{I1,I2};
    case 3
        found = false;
        for I1 = 1:length(fwd.table.grid)-1
            if fwd.table.grid(I1) <= m(1) && m(1) <= fwd.table.grid(I1+1)
                found = true;
                break
            end
        end
        assert(found,'m1 not found')

        found = false;
        for I2 = 1:length(fwd.table.grid)-1
            if fwd.table.grid(I2) <= m(2) && m(2) <= fwd.table.grid(I2+1)
                found = true;
                break
            end
        end
        assert(found,'m2 not found')

        found = false;
        for I3 = 1:length(fwd.table.grid)-1
            if fwd.table.grid(I3) <= m(3) && m(3) <= fwd.table.grid(I3+1)
                found = true;
                break
            end
        end
        assert(found,'m3 not found')

        obs = fwd.table.solVector{I1,I2,I3};
end
