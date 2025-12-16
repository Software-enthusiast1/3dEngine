# JS 3D Engine

This is a small 3D engine implemented entirely in `engine.js` and displayed by `index.html`.

**Controls**
- **Move:** `W` / `A` / `S` / `D`
- **Look (yaw/pitch):** Arrow keys
- **Roll:** `Q` (left) / `E` (right)
- **Jump:** Space (when grounded)
- **Reset camera:** `R`
- **Regenerate world (new seed):** `G`

**Run locally**
- Quick: open `index.html` directly in a modern browser (may have mixed results with modules).
- Recommended: run a local server (enables ES modules and CDN imports):

```bash
python3 -m http.server 8000
# then open http://localhost:8000
```

**What this engine does**
- Software renderer on an HTML canvas (triangle rasterization, back-face culling, simple lighting).
- Seeded procedural terrain using simplex-noise (fBm) with multiple octaves.
- Biomes: ocean, lake, desert, plains, snowy_plains, mountains (biome mix controlled by multi-scale noise).
- Vegetation: clustered tree placement with outliers and biome-aware types (cactus, oak, evergreen, shrub).

**Tuning & debugging**
- To reduce sharp ledges: increase smoothing iterations in `engine.js` (search `smoothIterations`) or lower noise frequency.
- To change water level / amplitude: edit the water scaling near `heightAt` in `engine.js`.