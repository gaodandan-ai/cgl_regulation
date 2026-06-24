import os
import re
import json
import math
import urllib.request
import urllib.error
import time

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
LITERATURE_DIR = os.path.join(BASE_DIR, "data", "literature")
CACHE_FILE = os.path.join(BASE_DIR, "data", "literature_cache.json")

class RAGService:
    def __init__(self):
        self.cache = {"files": {}}
        self.load_cache()
        self.ensure_dirs()

    def ensure_dirs(self):
        if not os.path.exists(LITERATURE_DIR):
            os.makedirs(LITERATURE_DIR)
        readme_path = os.path.join(LITERATURE_DIR, "README.txt")
        if not os.path.exists(readme_path):
            try:
                with open(readme_path, "w", encoding="utf-8") as f:
                    f.write("在此放置您的本地文献（支持 .txt, .md 格式）。\nRAG 知识库系统会自动扫描这些文件，提取关于基因或通路的调控规则，作为 AI 文献智能总结的补充背景知识。")
            except Exception as e:
                print("Error writing RAG README:", e)

    def load_cache(self):
        if os.path.exists(CACHE_FILE):
            try:
                with open(CACHE_FILE, "r", encoding="utf-8") as f:
                    self.cache = json.load(f)
                    if "files" not in self.cache:
                        self.cache = {"files": {}}
            except Exception as e:
                print("Error loading RAG cache:", e)
                self.cache = {"files": {}}

    def save_cache(self):
        try:
            with open(CACHE_FILE, "w", encoding="utf-8") as f:
                json.dump(self.cache, f, ensure_ascii=False, indent=2)
        except Exception as e:
            print("Error saving RAG cache:", e)

    def chunk_text(self, text, chunk_size=500, overlap=100):
        text = re.sub(r'\s+', ' ', text).strip()
        chunks = []
        start = 0
        while start < len(text):
            end = start + chunk_size
            if end >= len(text):
                chunks.append(text[start:])
                break
            
            boundary = -1
            for char in ['.', '!', '?', '。', '；', '\n']:
                pos = text.rfind(char, start + chunk_size - overlap, end)
                if pos > boundary:
                    boundary = pos
            
            if boundary != -1:
                end = boundary + 1
            
            chunks.append(text[start:end])
            start = end
        return chunks

    def fetch_embedding(self, text, provider, api_key, model_name, base_url):
        if provider == 'google':
            url = f"https://generativelanguage.googleapis.com/v1beta/models/text-embedding-004:embedContent?key={api_key}"
            payload = {
                "model": "models/text-embedding-004",
                "content": {"parts": [{"text": text}]}
            }
            try:
                req = urllib.request.Request(
                    url,
                    data=json.dumps(payload).encode('utf-8'),
                    headers={'Content-Type': 'application/json'},
                    method='POST'
                )
                with urllib.request.urlopen(req) as resp:
                    res = json.loads(resp.read().decode('utf-8'))
                    return res.get("embedding", {}).get("values", [])
            except Exception as e:
                print(f"Gemini embedding failed: {e}")
                return None
        elif provider == 'ollama':
            url = base_url if base_url else "http://localhost:11434"
            if not url.endswith('/api/embeddings') and not url.endswith('/v1/embeddings'):
                url = url.rstrip('/') + "/api/embeddings"
            
            payload = {
                "model": model_name if model_name else "deepseek-r1",
                "prompt": text
            }
            if "/v1/embeddings" in url:
                payload = {
                    "model": model_name if model_name else "deepseek-r1",
                    "input": text
                }
            try:
                req = urllib.request.Request(
                    url,
                    data=json.dumps(payload).encode('utf-8'),
                    headers={'Content-Type': 'application/json'},
                    method='POST'
                )
                with urllib.request.urlopen(req) as resp:
                    res = json.loads(resp.read().decode('utf-8'))
                    if "embedding" in res:
                        return res["embedding"]
                    elif "data" in res and len(res["data"]) > 0:
                        return res["data"][0].get("embedding", [])
            except Exception as e:
                print(f"Ollama embedding failed: {e}")
                return None
        else:
            url = base_url if base_url else "https://api.openai.com/v1"
            url = url.rstrip('/') + "/embeddings"
            headers = {
                'Content-Type': 'application/json',
            }
            if api_key:
                headers['Authorization'] = f"Bearer {api_key}"
            
            emb_model = "text-embedding-3-small"
            if provider == 'zhipu':
                emb_model = "embedding-2"
            elif provider == 'qwen':
                emb_model = "text-embedding-v2"

            payload = {
                "model": emb_model,
                "input": text
            }
            try:
                req = urllib.request.Request(
                    url,
                    data=json.dumps(payload).encode('utf-8'),
                    headers=headers,
                    method='POST'
                )
                with urllib.request.urlopen(req) as resp:
                    res = json.loads(resp.read().decode('utf-8'))
                    data = res.get("data", [])
                    if data:
                        return data[0].get("embedding", [])
            except Exception as e:
                print(f"OpenAI-compatible embedding failed ({provider}): {e}")
                return None

    def sync_literature(self, provider, api_key, model_name, base_url):
        if not os.path.exists(LITERATURE_DIR):
            return

        changed = False
        files_in_dir = [f for f in os.listdir(LITERATURE_DIR) if f.endswith(('.txt', '.md')) and f != "README.txt"]
        
        cached_files = list(self.cache["files"].keys())
        for f in cached_files:
            if f not in files_in_dir:
                del self.cache["files"][f]
                changed = True

        for filename in files_in_dir:
            file_path = os.path.join(LITERATURE_DIR, filename)
            mtime = os.path.getmtime(file_path)
            
            cached = self.cache["files"].get(filename)
            if cached and cached.get("mtime") == mtime:
                continue

            try:
                with open(file_path, "r", encoding="utf-8") as f:
                    text = f.read()
            except Exception as e:
                print(f"Error reading {filename}: {e}")
                continue

            chunks = self.chunk_text(text)
            chunk_objs = []
            
            print(f"Indexing RAG file {filename} ({len(chunks)} chunks)...")
            for i, chunk in enumerate(chunks):
                vector = None
                if api_key or provider == 'ollama':
                    vector = self.fetch_embedding(chunk, provider, api_key, model_name, base_url)
                
                chunk_objs.append({
                    "text": chunk,
                    "vector": vector
                })
            
            self.cache["files"][filename] = {
                "mtime": mtime,
                "chunks": chunk_objs
            }
            changed = True

        if changed:
            self.save_cache()

    def get_tfidf_similarity(self, query, chunk_text):
        def tokenize(text):
            words = re.findall(r'\b\w+\b', text.lower())
            return [w for w in words if len(w) > 1]
        
        q_words = tokenize(query)
        c_words = tokenize(chunk_text)
        
        if not q_words or not c_words:
            return 0.0

        q_counts = {}
        for w in q_words:
            q_counts[w] = q_counts.get(w, 0) + 1

        c_counts = {}
        for w in c_words:
            c_counts[w] = c_counts.get(w, 0) + 1

        dot_product = 0.0
        for w in q_counts:
            if w in c_counts:
                dot_product += q_counts[w] * c_counts[w]

        q_len = math.sqrt(sum(v**2 for v in q_counts.values()))
        c_len = math.sqrt(sum(v**2 for v in c_counts.values()))

        if q_len == 0 or c_len == 0:
            return 0.0
            
        return dot_product / (q_len * c_len)

    def query_similarity(self, query, provider, api_key, model_name, base_url, top_n=3):
        try:
            self.sync_literature(provider, api_key, model_name, base_url)
        except Exception as e:
            print("Error syncing literature in RAG query:", e)

        all_chunks = []
        for filename, file_data in self.cache["files"].items():
            for idx, chunk in enumerate(file_data.get("chunks", [])):
                all_chunks.append({
                    "file": filename,
                    "text": chunk["text"],
                    "vector": chunk.get("vector")
                })

        if not all_chunks:
            return []

        query_vector = None
        has_vectors = any(c["vector"] is not None for c in all_chunks)
        
        if has_vectors and (api_key or provider == 'ollama'):
            query_vector = self.fetch_embedding(query, provider, api_key, model_name, base_url)

        results = []
        if query_vector:
            for chunk in all_chunks:
                vec = chunk["vector"]
                if not vec or len(vec) != len(query_vector):
                    score = self.get_tfidf_similarity(query, chunk["text"])
                else:
                    dot = sum(a * b for a, b in zip(query_vector, vec))
                    a_len = math.sqrt(sum(a**2 for a in query_vector))
                    b_len = math.sqrt(sum(b**2 for b in vec))
                    score = dot / (a_len * b_len) if a_len > 0 and b_len > 0 else 0.0
                results.append((score, chunk))
        else:
            for chunk in all_chunks:
                score = self.get_tfidf_similarity(query, chunk["text"])
                results.append((score, chunk))

        results.sort(key=lambda x: x[0], reverse=True)
        top_results = []
        for score, chunk in results[:top_n]:
            if score > 0.05:
                top_results.append({
                    "file": chunk["file"],
                    "text": chunk["text"],
                    "score": score
                })
        return top_results
