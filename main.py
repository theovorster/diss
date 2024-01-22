from flask import Flask, render_template, redirect, request, jsonify, send_file
from flask_bootstrap import Bootstrap5
from flask_wtf import FlaskForm
from wtforms import SubmitField, IntegerField, DecimalField, RadioField
from wtforms.validators import DataRequired, NumberRange
import os
import shutil
import pandas as pd


# Own python files in directory
from initial_processing import initial_processing
from secondary_processing import secondary_processing, tertiary_processing


UPLOAD_FOLDER = 'uploads'
iupred_number = None


app = Flask(__name__)
app.config['SECRET_KEY'] = '8BYkEfBA6O6donzWlSihBXox7C0sKR6b'
Bootstrap5(app)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


class ParameterForm(FlaskForm):
    disorder_score = DecimalField("Min. Disorder Score (0 - 1)", validators=[DataRequired(), NumberRange(min=0, max=1, message="Please give a number between 0 and 1.")])
    sequence_length = IntegerField("Min. Sequence Length (1 - ∞)", validators=[DataRequired()])
    merge_closer_than = IntegerField("Merge Sequences if Closer than X (1 - ∞)", validators=[DataRequired()])
    nuclear_score = DecimalField("Min. Nuclear Score (0 - 1)", validators=[DataRequired(), NumberRange(min=0, max=1, message="Please give a number between 0 and 1.")])
    iupred_or_anchor = RadioField("Search by IUPred or Anchor Score?", choices=["IUPred", "Anchor"], validators=[DataRequired()])
    submit = SubmitField('Submit')


@app.route('/', methods=['GET', 'POST'])
def home():
    return render_template("index.html")


@app.route('/how-to-use/', methods=['GET','POST'])
def how_to_use():
    return render_template("use_instructions.html")


@app.route('/upload')
def upload_page():
    return render_template('upload.html')


@app.route('/process_upload', methods=['POST'])
def process_upload():
    global iupred_number
    iupred_results = request.files['iupredResults']
    nuclear_scores = request.files['nuclearScores']

    if iupred_results.filename == '' or nuclear_scores.filename == '':
        return jsonify({'error': 'Please select both files'}), 400

    # Delete existing files
    try:
       shutil.rmtree(UPLOAD_FOLDER)
    except FileNotFoundError:
        pass
    os.mkdir(UPLOAD_FOLDER)

    # Save new files
    iupred_results.save(os.path.join(app.config['UPLOAD_FOLDER'], 'iupred.txt'))
    nuclear_scores.save(os.path.join(app.config['UPLOAD_FOLDER'], 'nuclear.csv'))

    # Process iupred.txt and nuclear.csv
    iupred_number = initial_processing(UPLOAD_FOLDER)
    

    return redirect("/process-files")


@app.route("/process-files", methods=["GET", "POST"])
def process_files():
    form = ParameterForm()

    # Form for processing parameters.
    if form.validate_on_submit():
        disorder = float(form.disorder_score.data)
        sequence = int(form.sequence_length.data)
        merge = int(form.merge_closer_than.data)
        nuclear_score = float(form.nuclear_score.data)
        iupred_or_anchor = form.iupred_or_anchor.data
        print(nuclear_score)
        print(iupred_or_anchor)
        print("wah")
        secondary_processing(iupred_number, disorder, sequence, merge, iupred_or_anchor)
        tertiary_processing(nuclear_score)
        return redirect("/long-dataset")

    return render_template("process_files.html", form=ParameterForm())


@app.route("/long-dataset", methods=["GET", "POST"])
def dataset():
    # Reads and displays the completed table.

    try:
        df = pd.read_csv(
            "data/final_results/long_dataset.csv",
            encoding="unicode-escape",
            usecols=[
                "Identifier",
                "Disordered region",
                "Disorder score",
                "SIM Position Site",
                "SIM Sequence",
                "SIM Type",
                "SIM region sequence",
                "D/E",
                "S/T",
                "P",
                "Nuclear Score"
                ]
            )
        df = df.reindex(columns=[
                "Identifier",
                "Disordered region",
                "Disorder score",
                "Nuclear Score",
                "SIM Position Site",
                "SIM Sequence",
                "SIM Type",
                "SIM region sequence",
                "D/E",
                "S/T",
                "P"
                ])

        # Create a list of dictionaries to pass to the template
        data_list = df.to_dict(orient='records')
        
        return render_template("long_dataset.html", data=data_list, csv_found=True)
    except FileNotFoundError:
        return render_template("long_dataset.html", csv_found=False)


@app.route("/download_long_csv", methods=["POST"])
def download_long_csv():
    file_path = "data/final_results/long_dataset.csv"
    return send_file(file_path, as_attachment=True)


@app.route("/short_dataset", methods=["GET", "POST"])
def short_dataset():
    # Reads and displays the completed table.
    try:
        df = pd.read_csv(
            "data/final_results/short_dataset.csv",
            encoding="unicode-escape",
            )
        df = df.sort_values(by=["No. Putative SIMs"], ascending=False)
        return render_template("short_dataset.html", data=df.to_html(classes="table table-hover", index=False, justify="center"), csv_found=True)
    except FileNotFoundError:
        return render_template("short_dataset.html", csv_found=False)


@app.route("/download_short_csv", methods=["POST"])
def download_short_csv():
    file_path = "data/final_results/short_dataset.csv"
    return send_file(file_path, as_attachment=True)


@app.route("/about")
def about():
    return render_template("about.html")


@app.route("/contact")
def contact():
    return render_template("contact.html")


if __name__ == "__main__":
    app.run(debug=True, port=5002)
