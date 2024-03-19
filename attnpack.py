import os
from protein_learning.models.inference_utils import Inference, make_predicted_protein
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputpdb", required=True)
    parser.add_argument("--outputpath", required=True)
    parser.add_argument("--numcores")
    args = vars(parser.parse_args())
    numcores = args["numcores"] if args["numcores"] else "4"
    os.environ["OPENBLAS_NUM_THREADS"] = numcores
    os.environ["MKL_NUM_THREADS"] = numcores
    os.environ["OMP_NUM_THREADS"] = numcores

    RESOURCE_ROOT = "/data/jgut/template-analysis/AttnPacker/AttnPackerPTM_V2"
    runner = Inference(RESOURCE_ROOT, use_design_variant=False).to("cpu")
    prediction = runner.infer(
        pdb_path=args["inputpdb"],
        # Boolean Tensor indicating which residues to design
        seq_mask=None,
        # Whether to format output (process into logits, seq labels, etc. or return raw node and pair output)
        format=True,
        # Chunk inference by successively packing smaller crops of size chunk_size
        # This allows packing of arbitrarily long proteins
        chunk_size=500,
    )

    predicted_protein = make_predicted_protein(
        model_out=prediction["model_out"], seq=prediction["seq"]
    )
    # write predicted protein with plddt score in b-factor column
    predicted_protein.to_pdb(
        args["outputpath"],
        beta=[elem * 100 for elem in prediction["pred_plddt"].squeeze()],
    )
