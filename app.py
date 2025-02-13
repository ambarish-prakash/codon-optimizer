from flask import Flask, request, render_template
from dnachisel import *
import pandas as pd
import os

app = Flask(__name__)
PORT = int(os.getenv("PORT", 5000)) 

def load_codon_usage_table(csv_path):
    df = pd.read_csv(csv_path)
    df["codon"] = df["codon"].str.replace("U", "T")  # Convert RNA codons (U → T)

    codon_table = (
        df.groupby("amino_acid", group_keys=False)
        .apply(lambda x: dict(zip(x["codon"], x["relative_frequency"])))
        .to_dict()
    )
    
    return codon_table

codon_usage_table = load_codon_usage_table("k_rhaeticus.csv")

def optimize_sequence(sequence):
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[
            AvoidPattern('BbsI_site'),
            AvoidPattern('BsmBI_site'),
            AvoidPattern('BsaI_site'),
            AvoidPattern('NotI_site'),
            EnforceTranslation(genetic_table="Bacterial"),
            EnforceGCContent(mini=0.3, maxi=0.7, window=50),
        ],
        objectives=[CodonOptimize(codon_usage_table=codon_usage_table)]
    )

    problem.resolve_constraints()
    problem.optimize()

    return {
        "optimized_sequence": problem.sequence,
        "constraints_summary": problem.constraints_text_summary(),
        "objectives_summary": problem.objectives_text_summary()
    }

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        sequence = request.form['sequence']
        results = optimize_sequence(sequence)

        return render_template('index.html', 
            input_sequence=sequence, 
            optimized_sequence=results["optimized_sequence"],
            constraints_summary=results["constraints_summary"],
            objectives_summary=results["objectives_summary"]
        )
    
    return render_template('index.html', input_sequence="", optimized_sequence="", constraints_summary="", objectives_summary="")

if __name__ == '__main__':
    app.run(host="0.0.0.0", port=PORT)
