function tile = combineTiles(tile0, tile180)

tile = tile0;

tile.coords = [tile0.coords; tile180.coords];
tile.taps   = [tile0.taps;   tile180.taps];
tile.mean   = [tile0.mean;   tile180.mean];
tile.std    = [tile0.std;    tile180.std];
tile.peak   = [tile0.peak';   tile180.peak'];

if tile0.CI == 'on'
    tile.CI95_mean = [tile0.CI95_mean; tile180.CI95_mean];
    tile.CI95_std  = [tile0.CI95_std;  tile180.CI95_std];
    tile.CI95_peak = [tile0.CI95_peak;  tile180.CI95_peak];
end
