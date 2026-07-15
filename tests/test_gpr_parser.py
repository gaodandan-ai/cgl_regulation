import sys
import os

# Add backend to sys.path so we can import gene_utils
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "backend"))

from gene_utils import evaluate_gpr_rule

def test_single_gene():
    gene_values = {"Cgl0851": 2.5}
    assert evaluate_gpr_rule("Cgl0851", gene_values) == 2.5
    assert evaluate_gpr_rule("Cgl0000", gene_values, default_val=1.0) == 1.0

def test_and_operator():
    gene_values = {"g1": 0.5, "g2": 2.0}
    assert evaluate_gpr_rule("g1 and g2", gene_values) == 0.5
    assert evaluate_gpr_rule("g1 AND g2", gene_values) == 0.5

def test_or_operator():
    gene_values = {"g1": 0.5, "g2": 2.0}
    assert evaluate_gpr_rule("g1 or g2", gene_values) == 2.0
    assert evaluate_gpr_rule("g1 OR g2", gene_values) == 2.0

def test_nested_expression():
    gene_values = {"g1": 0.5, "g2": 2.0, "g3": 1.5}
    # (g1 or g2) and g3 -> max(0.5, 2.0) and 1.5 -> 2.0 and 1.5 -> min(2.0, 1.5) = 1.5
    assert evaluate_gpr_rule("(g1 or g2) and g3", gene_values) == 1.5
    
    # g1 or (g2 and g3) -> 0.5 or min(2.0, 1.5) -> 0.5 or 1.5 -> max(0.5, 1.5) = 1.5
    assert evaluate_gpr_rule("g1 or (g2 and g3)", gene_values) == 1.5

def test_missing_gene_fallback():
    gene_values = {"g1": 0.5}
    # g1 and g_missing -> min(0.5, 1.0) = 0.5
    assert evaluate_gpr_rule("g1 and g_missing", gene_values, default_val=1.0) == 0.5
    # g1 or g_missing -> max(0.5, 1.0) = 1.0
    assert evaluate_gpr_rule("g1 or g_missing", gene_values, default_val=1.0) == 1.0

if __name__ == "__main__":
    test_single_gene()
    test_and_operator()
    test_or_operator()
    test_nested_expression()
    test_missing_gene_fallback()
    print("All evaluate_gpr_rule tests passed successfully!")
