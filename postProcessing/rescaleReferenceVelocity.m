function tile = rescaleReferenceVelocity(tile_old, Uref)

tile = tile_old;

tile.mean = tile_old.mean * (0.5*tile.U_1m^2) / (0.5*Uref^2);
tile.std  = tile_old.std  * (0.5*tile.U_1m^2) / (0.5*Uref^2);
tile.peak = tile_old.peak * (0.5*tile.U_1m^2) / (0.5*Uref^2);

% tile.timeHistory = tile_old.timeHistory * (0.5*tile.U_1m^2) / (0.5*Uref^2);
% tile.areaAverage = tile_old.areaAverage * (0.5*tile.U_1m^2) / (0.5*Uref^2);

tile.U = Uref;





