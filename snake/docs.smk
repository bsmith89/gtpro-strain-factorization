rule render_figure_to_png:
    output: "fig/{stem}_figure.w{width}.png"
    input: "doc/static/{stem}_figure.svg"
    params:
        width=lambda w: int(w.width)
    shell:
        """
        inkscape {input} --export-width={params.width} --export-filename {output}
        """
