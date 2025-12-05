# JS 3D Engine With Terrain Gen.

This is a 3D engine implemented entirely inside `engine.js` and displayed by `index.html`.

Controls:
- `W` / `A` / `S` / `D` — move camera forward/left/back/right
- Arrow keys — rotate camera (yaw & pitch)
- `R` — reset camera position
- `G` — generate new world

To run:
1. Open `index.html` in your browser (double-click or run a local HTTP server).
2. For best results run a simple local server from the project folder, e.g.:

```bash
python3 -m http.server 8000
# then open http://localhost:8000 in your browser
```

The engine uses a software renderer on a 2D canvas (filled triangles, back-face culling and simple lighting) made from scratch. The seed based terrain generator uses perlin noise, is also made from scratch

Please tell me how I can fix the terrain generaton because it likes to sometimes generate sharp ledges.