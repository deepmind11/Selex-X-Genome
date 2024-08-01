from argparse import ArgumentParser
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from ENCODE.experiment import TFChipSeq
from matplotlib.backends.backend_pdf import PdfPages


def list_png_files(directory):
    try:
        directory_path = Path(directory)
        if not directory_path.exists():
            return f"The directory {directory} does not exist."
        png_files = [
            str(file) for file in directory_path.glob("*FOXA1_X_FOXA1_Enr_vs_Bin.png")
        ]
        return png_files
    except Exception as e:
        return str(e)


if __name__ == "__main__":

    parser = ArgumentParser()

    # Required positional arguments
    parser.add_argument("TF_Directory", help="Path to all the experiments")
    args = parser.parse_args()

    TF_Directory = Path(args.TF_Directory)
    # Path to the pdf
    pdf_foxa1 = Path("/burg/home/hg2604/FOXA1_plot.pdf")

    # All the experiments
    experiment_drs = [entry for entry in TF_Directory.iterdir() if entry.is_dir()]

    with PdfPages(pdf_foxa1) as pdf:

        for experiment in experiment_drs:

            # Insert a blank page at the start of each experiment.
            expr_acc = experiment.name
            encode_experiment_obj = TFChipSeq(expr_acc)
            meta_data = encode_experiment_obj.get_other_meta_data()

            description = f"Summary: {encode_experiment_obj.expr_data['simple_biosample_summary']}"
            description += "\n"
            for key in meta_data.keys():
                description += f"{key}: {meta_data[key]}"
                description += "\n"

            library_drs = [entry for entry in experiment.iterdir() if entry.is_dir()]

            for library in library_drs:

                png_files = list_png_files(library)

                for image in png_files:

                    title = f"{expr_acc}-{library.name}-{Path(image).name}"

                    img = mpimg.imread(image)
                    plt.imshow(img)
                    plt.title(title)
                    plt.axis("off")  # Optional: turn off the axis

                    # Add description at the bottom
                    plt.figtext(
                        0.5, 0.02, description, ha="center", fontsize=10, wrap=True
                    )

                    # Adjust the layout to make space for the description
                    plt.subplots_adjust(bottom=0.5)

                    pdf.savefig()  # saves the current figure into a pdf page
                    plt.close()
