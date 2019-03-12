"""Write report about the model run to a Jupyter Notebook

Description:
------------

asdf

"""
# Standard library imports
import itertools

# WHERE imports
from where.lib import files
from where.lib import log
from where.lib import plugins
from where.reports import report

# Optional library imports
from midgard.dev import optional

kernelspec = optional.optional_import("jupyter_client.kernelspec")
nbformat = optional.optional_import("nbformat")
nbbase = optional.optional_import("nbformat.v4.nbbase")

CELL_PARSER = {
    "title": lambda txt: "# " + txt,
    "header": lambda txt: "## " + txt,
    "text": lambda txt: txt,
    "model": lambda txt: "### Model: " + txt.replace("_", " ").title(),
}


@plugins.register
def write_notebook(rundate, tech):
    kernelspec = _get_kernelspec("python3")
    cells = list(get_cells())

    notebook = nbbase.new_notebook(cells=cells, metadata={"language": "python", "kernelspec": kernelspec})

    nb_filename = files.path("output_notebook")
    nbformat.write(notebook, nb_filename)
    log.info(f"Create Jupyter Notebook report. To open:\n        jupyter notebook {nb_filename}")


def _get_kernelspec(name):
    ksm = kernelspec.KernelSpecManager()
    ks_dict = ksm.get_kernel_spec(name).to_dict()
    ks_dict["name"] = name
    ks_dict.pop("argv")

    return ks_dict


def get_cells():
    slide_type = itertools.chain(["slide", "fragment"], itertools.cycle(["subslide", "fragment"]))
    for report_name, _, report_data in report.reports():
        #        print(report_name, report_data)
        _, report_type = report_name.split("/")
        if report_type == "header":
            slide_type = itertools.chain(["slide", "fragment"], itertools.cycle(["subslide", "fragment"]))
        if report_type == "model":
            slide_type = itertools.cycle(["subslide", "fragment"])
        if "__doc__" in report_data and report_data["__doc__"]:
            yield nbbase.new_markdown_cell(report_data["__doc__"], metadata=_meta_slideshow("fragment"))
        if "__code__" in report_data:
            yield nbbase.new_code_cell(report_data["__code__"], metadata=_meta_slideshow("fragment"))
        else:
            parser = CELL_PARSER.get(report_type)
            if parser:
                yield nbbase.new_markdown_cell(parser(report_data["text"]), metadata=_meta_slideshow(next(slide_type)))


def _meta_slideshow(slide_type):
    return dict(slideshow={"slide_type": slide_type})
