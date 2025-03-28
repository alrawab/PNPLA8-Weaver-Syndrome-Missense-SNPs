# ==============================================================================================
#                  invoke_missense_SNPs_PNPLA8_Weaver_Syndrome.py
# 
# Description: 
# This Python script is designed to identify and analyze missense SNPs in the PNPLA8 gene 
# that are associated with Weaver Syndrome in cattle. The script processes genomic data, 
# identifies mutations, and evaluates their potential impact on protein function.
#
# Author: Dr. Osamah S. Alraouwab
# Department: Molecular Biology and Biochemistry
# Institution: Zintan Faculty of Medicine
#
# Contact Information:
# Email: rawab@uoz.edu.ly
# 
#
# Notes:
# - Ensure that all required libraries and dependencies are installed before running the script.
# - Input data should be in a compatible format (e.g., VCF, FASTA) for proper analysis.
# - The output will include detailed information about identified missense SNPs and their 
#   potential effects on the PNPLA8 gene.
#
# Disclaimer:
# This script is intended for research purposes only. Results should be validated through 
# additional experimental methods.
#
# License:
# This program is free software: you can redistribute it and/or modify it under the terms of 
# the GNU General Public License as published by the Free Software Foundation, either version 3 
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. 
# If not, see <https://www.gnu.org/licenses/>.
# ==============================================================================================
import requests
import pandas as pd
from pathlib import Path
import time
import json
import numpy as np
from bs4 import BeautifulSoup
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configuration
output_dir = Path("db")
output_dir.mkdir(exist_ok=True)
output_csv = output_dir / "pnpla8_missense_enhanced.csv"

chromosome = "4"
start = 49_588_563
end = 49_654_963
API_DELAY = 0.5  # Seconds between API calls
MAX_RETRIES = 3
MAX_WORKERS = 4  # For parallel processing

def classify_significance(row):
    """Enhanced clinical significance classification"""
    sift = str(row['sift_pred']).lower() if pd.notna(row['sift_pred']) else ''
    polyphen = str(row['polyphen_pred']).lower() if pd.notna(row['polyphen_pred']) else ''
    
    sift_deleterious = 'deleterious' in sift
    polyphen_damaging = any(x in polyphen for x in ['damaging', 'probably'])
    sift_tolerated = 'tolerated' in sift
    polyphen_benign = 'benign' in polyphen
    
    if sift_deleterious and polyphen_damaging:
        return "Pathogenic"
    elif sift_deleterious or polyphen_damaging:
        return "Likely Pathogenic"
    elif sift_tolerated and polyphen_benign:
        return "Likely Benign"
    elif sift_tolerated or polyphen_benign:
        return "Possibly Benign"
    return "VUS"

def get_polyphen_prediction(aa_change):
    """Get PolyPhen-2 predictions via web API with retries"""
    if not aa_change or pd.isna(aa_change) or not isinstance(aa_change, str):
        return None, None
    
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.post(
                "http://genetics.bwh.harvard.edu/pph2/bgi.cgi",
                data={'sequence': f">variant\n{aa_change}", 'submit': 'Submit'},
                timeout=15
            )
            soup = BeautifulSoup(response.text, 'html.parser')
            result = soup.find('pre').text
            
            if "PolyPhen-2 score:" not in result:
                raise ValueError("PolyPhen response format unexpected")
                
            score = float(result.split('PolyPhen-2 score:')[1].split()[0])
            pred = result.split('Prediction:')[1].split()[0].lower()
            return score, pred
            
        except Exception as e:
            print(f"PolyPhen attempt {attempt + 1} for {aa_change} failed: {str(e)}")
            time.sleep(API_DELAY * (attempt + 1))
    return None, None

def parse_amino_acid_change(hgvsp):
    """Extract reference and alternate amino acids from HGVSp notation"""
    if not hgvsp or not isinstance(hgvsp, str):
        return None, None, None
    
    try:
        # Example: p.Arg345Cys → ('Arg', '345', 'Cys')
        parts = hgvsp.split(':')[-1].replace('p.', '').split('.')[-1]
        ref_aa = ''.join([c for c in parts if c.isalpha() or c == '_']).split('_')[0]
        pos = ''.join([c for c in parts if c.isdigit()])
        alt_aa = ''.join([c for c in parts[::-1] if c.isalpha() or c == '_'])[::-1].split('_')[0]
        
        return ref_aa, pos, alt_aa
    except:
        return None, None, None

def get_vep_annotation(variant_id):
    """Get variant annotations from Ensembl VEP with enhanced PolyPhen handling"""
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(
                f"https://rest.ensembl.org/vep/bos_taurus/id/{variant_id}",
                headers={"Content-Type": "application/json"},
                timeout=30
            )
            
            if not response.text.strip():
                print(f"Empty response for {variant_id} (attempt {attempt + 1})")
                time.sleep(API_DELAY * (attempt + 1))
                continue
                
            data = response.json()
            if not data:
                return None
                
            for item in data:
                if "transcript_consequences" in item:
                    for tc in item["transcript_consequences"]:
                        if "sift_prediction" in tc:
                            # Parse amino acid change
                            hgvsp = tc.get("hgvsp", "")
                            ref_aa, aa_pos, alt_aa = parse_amino_acid_change(hgvsp)
                            
                            # Get standard annotations
                            annotation = {
                                "amino_acid_change": hgvsp.split(':')[-1] if hgvsp else None,
                                "ref_aa": ref_aa,
                                "alt_aa": alt_aa,
                                "aa_position": aa_pos,
                                "protein_position": tc.get("protein_start"),
                                "alleles": f"{tc.get('allele_string','')}",
                                "consequence": ", ".join(tc.get("consequence_terms", [])),
                                "sift_score": tc.get("sift_score"),
                                "sift_pred": tc.get("sift_prediction"),
                                "polyphen_score": None,
                                "polyphen_pred": None,
                                "polyphen_source": "None",
                                "confidence": "low" if "low_confidence" in str(tc.get("sift_prediction", "")).lower() else "high"
                            }
                            
                            # Try to get PolyPhen data from multiple possible fields
                            polyphen = tc.get("polyphen", {})
                            if isinstance(polyphen, dict):
                                annotation.update({
                                    "polyphen_score": polyphen.get("score"),
                                    "polyphen_pred": polyphen.get("prediction"),
                                    "polyphen_source": "Ensembl" if polyphen.get("score") else "None"
                                })
                            else:
                                annotation.update({
                                    "polyphen_score": tc.get("polyphen_score"),
                                    "polyphen_pred": tc.get("polyphen_prediction"),
                                    "polyphen_source": "Ensembl" if tc.get("polyphen_score") else "None"
                                })
                            
                            # If still missing, try external PolyPhen-2
                            if pd.isna(annotation["polyphen_pred"]):
                                aa_change = annotation["amino_acid_change"]
                                if aa_change and aa_change != "":
                                    score, pred = get_polyphen_prediction(aa_change)
                                    annotation.update({
                                        "polyphen_score": score,
                                        "polyphen_pred": pred,
                                        "polyphen_source": "PolyPhen-2" if score else "None"
                                    })
                            
                            return annotation
            return None
            
        except json.JSONDecodeError:
            print(f"Invalid JSON for {variant_id} (attempt {attempt + 1})")
            time.sleep(API_DELAY * (attempt + 1))
        except Exception as e:
            print(f"Error processing {variant_id} (attempt {attempt + 1}): {str(e)}")
            time.sleep(API_DELAY * (attempt + 1))
    return None

def process_variant(var):
    """Process a single variant with error handling"""
    variant_id = var["id"]
    print(f"Processing {variant_id}")
    
    annotation = get_vep_annotation(variant_id)
    if annotation:
        result = {
            "rs_id": variant_id,
            "chromosome": chromosome,
            "position": var["start"],
            **annotation
        }
        time.sleep(API_DELAY)
        return result
    return None

# Main execution
if __name__ == "__main__":
    print("Fetching variants from Ensembl...")
    response = requests.get(
        f"https://rest.ensembl.org/overlap/region/Bos_taurus/{chromosome}:{start}-{end}"
        "?feature=variation;content-type=application/json"
    )
    all_variants = response.json()

    missense_variants = [
        var for var in all_variants 
        if "missense_variant" in var.get("consequence_type", [])
    ]
    print(f"Found {len(missense_variants)} potential missense variants")

    # Process variants in parallel
    results = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = [executor.submit(process_variant, var) for var in missense_variants]
        for future in as_completed(futures):
            result = future.result()
            if result:
                results.append(result)

    # Create DataFrame and enhance
    if results:
        df = pd.DataFrame(results)
        
        # Extract protein position if missing
        if 'protein_position' not in df.columns:
            df['protein_position'] = df['amino_acid_change'].str.extract(r'(\d+)')
        
        # Add clinical significance
        df['clinical_significance'] = df.apply(classify_significance, axis=1)
        
        # Final column order
        column_order = [
            'rs_id', 'chromosome', 'position', 
            'alleles', 'ref_aa', 'alt_aa', 'aa_position',
            'amino_acid_change', 'protein_position',
            'consequence', 'clinical_significance',
            'sift_score', 'sift_pred', 
            'polyphen_score', 'polyphen_pred', 'polyphen_source',
            'confidence'
        ]
        df = df[column_order]
        
        # Save to CSV
        df.to_csv(output_csv, index=False)
        print(f"Successfully exported {len(df)} variants to {output_csv}")
        
        # Verification output
        print("\nFirst 5 variants:")
        print(df.head().to_string(index=False))
        
        # Summary stats
        print("\nClinical Significance Summary:")
        print(df['clinical_significance'].value_counts())
        
        print("\nPolyPhen Prediction Sources:")
        print(df['polyphen_source'].value_counts())
    else:
        print("No variants were successfully annotated")