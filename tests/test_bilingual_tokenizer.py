import sys
import os
import re

# Add root folder to sys.path so we can import rag_service
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rag_service import RAGService

def test_english_tokenization():
    rag = RAGService()
    # Simple tfidf logic helper
    sim_func = rag.get_tfidf_similarity
    
    # Let's test the nested tokenize function inside get_tfidf_similarity
    # We can retrieve it by parsing a query and seeing which terms overlap
    # Or we can just mock/test get_tfidf_similarity outputs.
    
    # If we query 'glutamate', it should score > 0 on English documents containing 'glutamate'
    score_match = sim_func("glutamate pathway", "This is the glutamate biosynthesis pathway.")
    score_mismatch = sim_func("glutamate pathway", "Lysine biosynthesis is highly regulated.")
    assert score_match > 0.1
    assert score_mismatch == 0.0

def test_chinese_tokenization():
    rag = RAGService()
    sim_func = rag.get_tfidf_similarity
    
    # Query '谷氨酸' (glutamate) in Chinese
    # If the tokenizer splits properly into characters ('谷', '氨', '酸') and bigrams ('谷氨', '氨酸'),
    # the overlap score with a document containing '谷氨酸合成' should be high.
    score_match = sim_func("谷氨酸", "这是谷氨酸合成通路。")
    assert score_match > 0.2
    
    # Query '赖氨酸' (lysine)
    score_diff = sim_func("赖氨酸", "这是谷氨酸合成通路。")
    # '赖氨酸' -> ('赖', '氨', '酸', '赖氨', '氨酸')
    # '谷氨酸' -> ('谷', '氨', '酸', '谷氨', '氨酸')
    # They overlap on '氨', '酸', '氨酸', so similarity score will be positive but lower than a perfect match.
    assert score_diff > 0.0
    assert score_match > score_diff

def test_mixed_bilingual_tokenization():
    rag = RAGService()
    sim_func = rag.get_tfidf_similarity
    
    # English word mixed with Chinese characters
    score_mix = sim_func("cg0251 基因", "这个 locus 是 cg0251 编码的蛋白。")
    assert score_mix > 0.1

if __name__ == "__main__":
    test_english_tokenization()
    test_chinese_tokenization()
    test_mixed_bilingual_tokenization()
    print("All bilingual tokenizer tests passed successfully!")
