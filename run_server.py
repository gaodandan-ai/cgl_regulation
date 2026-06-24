#!/usr/bin/env python3
import os
import http.server
import socketserver
import webbrowser
import threading
import time
import sys
import urllib.request
import urllib.parse
import json
import re
import urllib.error

PORT = 8000

# Caches for KEGG pathways and GO terms
KEGG_PATHWAY_NAMES = {}       # cgb/cgl pathway ID -> clean name
PATHWAY_NAMES_MUTEX = threading.Lock()
GENE_PATHWAYS_CACHE = {}      # (cg_locus, cgl_locus) -> parsed dict

def load_kegg_pathway_names():
    global KEGG_PATHWAY_NAMES
    if KEGG_PATHWAY_NAMES:
        return
        
    with PATHWAY_NAMES_MUTEX:
        if KEGG_PATHWAY_NAMES:
            return
            
        # Fetch Bielefeld (cgb) pathways
        try:
            url_cgb = "https://rest.kegg.jp/list/pathway/cgb"
            req_cgb = urllib.request.Request(url_cgb, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req_cgb, timeout=10) as resp:
                lines = resp.read().decode('utf-8').splitlines()
                for line in lines:
                    if '\t' in line:
                        pid, pname = line.split('\t', 1)
                        pname_clean = pname.split(" - Corynebacterium")[0].strip()
                        pid_clean = pid.strip()
                        KEGG_PATHWAY_NAMES[pid_clean] = pname_clean
                        if not pid_clean.startswith("path:"):
                            KEGG_PATHWAY_NAMES[f"path:{pid_clean}"] = pname_clean
        except Exception as e:
            print("Error loading cgb pathway names:", e)
            
        # Fetch Kyowa Hakko (cgl) pathways
        try:
            url_cgl = "https://rest.kegg.jp/list/pathway/cgl"
            req_cgl = urllib.request.Request(url_cgl, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req_cgl, timeout=10) as resp:
                lines = resp.read().decode('utf-8').splitlines()
                for line in lines:
                    if '\t' in line:
                        pid, pname = line.split('\t', 1)
                        pname_clean = pname.split(" - Corynebacterium")[0].strip()
                        pid_clean = pid.strip()
                        KEGG_PATHWAY_NAMES[pid_clean] = pname_clean
                        if not pid_clean.startswith("path:"):
                            KEGG_PATHWAY_NAMES[f"path:{pid_clean}"] = pname_clean
        except Exception as e:
            print("Error loading cgl pathway names:", e)

def get_gene_pathways_and_go(cg, cgl):
    global GENE_PATHWAYS_CACHE
    
    # Standardize tags
    cg = cg.strip() if cg else ""
    cgl = cgl.strip() if cgl else ""
    
    cache_key = (cg.lower(), cgl.lower())
    if cache_key in GENE_PATHWAYS_CACHE:
        return GENE_PATHWAYS_CACHE[cache_key]
        
    load_kegg_pathway_names()
    
    pathways = []
    seen_pids = set() # Store numeric part of pathway IDs (e.g. "02020")
    
    # 1. Query cgb pathways for cg_locus
    if cg:
        try:
            url = f"https://rest.kegg.jp/link/pathway/cgb:{cg}"
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=10) as resp:
                lines = resp.read().decode('utf-8').splitlines()
                for line in lines:
                    if '\t' in line:
                        _, pid_raw = line.split('\t', 1)
                        pid_clean = pid_raw.replace("path:", "").strip()
                        pid_num = "".join(c for c in pid_clean if c.isdigit())
                        if pid_num not in seen_pids:
                            seen_pids.add(pid_num)
                            name = KEGG_PATHWAY_NAMES.get(pid_clean, pid_clean)
                            link = f"https://www.kegg.jp/kegg-bin/show_pathway?{pid_clean}+cgb:{cg}"
                            pathways.append({
                                "id": pid_clean,
                                "name": name,
                                "link": link,
                                "source": "KEGG"
                            })
        except Exception as e:
            print(f"Error querying cgb pathways for {cg}: {e}")
            
    # 2. Query cgl pathways for cgl_locus
    if cgl:
        # Standardize capitalization: e.g. cgl0339 -> Cgl0339
        cgl_normalized = cgl
        if len(cgl) > 3 and cgl.lower().startswith('cgl'):
            cgl_normalized = 'Cgl' + cgl[3:]
            
        try:
            url = f"https://rest.kegg.jp/link/pathway/cgl:{cgl_normalized}"
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=10) as resp:
                lines = resp.read().decode('utf-8').splitlines()
                for line in lines:
                    if '\t' in line:
                        _, pid_raw = line.split('\t', 1)
                        pid_clean = pid_raw.replace("path:", "").strip()
                        pid_num = "".join(c for c in pid_clean if c.isdigit())
                        if pid_num not in seen_pids:
                            seen_pids.add(pid_num)
                            name = KEGG_PATHWAY_NAMES.get(pid_clean, pid_clean)
                            link = f"https://www.kegg.jp/kegg-bin/show_pathway?{pid_clean}+cgl:{cgl_normalized}"
                            pathways.append({
                                "id": pid_clean,
                                "name": name,
                                "link": link,
                                "source": "KEGG"
                            })
        except Exception as e:
            print(f"Error querying cgl pathways for {cgl_normalized}: {e}")
            
    # 3. Query GO terms from UniProt
    go_terms = []
    seen_gos = set()
    
    # Query UniProt using cg or cgl locus
    query_tag = cg if cg else cgl
    if query_tag:
        try:
            uniprot_url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{query_tag}+AND+organism_id:196627&fields=id,accession,go&format=json"
            req = urllib.request.Request(uniprot_url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=10) as resp:
                data = json.loads(resp.read().decode('utf-8'))
                results = data.get("results", [])
                if results:
                    refs = results[0].get("uniProtKBCrossReferences", [])
                    for ref in refs:
                        if ref.get("database") == "GO":
                            go_id = ref.get("id")
                            props = ref.get("properties", [])
                            go_term_val = ""
                            for prop in props:
                                if prop.get("key") == "GoTerm":
                                    go_term_val = prop.get("value")
                                    break
                            if go_id and go_term_val and go_id not in seen_gos:
                                seen_gos.add(go_id)
                                go_type = "GO"
                                go_name = go_term_val
                                if ":" in go_term_val:
                                    t_code, t_name = go_term_val.split(":", 1)
                                    if t_code == "P":
                                        go_type = "GO Process"
                                    elif t_code == "F":
                                        go_type = "GO Function"
                                    elif t_code == "C":
                                        go_type = "GO Component"
                                    go_name = t_name.strip()
                                
                                link = f"https://www.ebi.ac.uk/QuickGO/term/{go_id}"
                                go_terms.append({
                                    "id": go_id,
                                    "name": go_name,
                                    "type": go_type,
                                    "link": link
                                })
        except Exception as e:
            print(f"Error querying UniProt GO terms for {query_tag}: {e}")
            
    result = {
        "pathways": pathways,
        "go_terms": go_terms
    }
    
    GENE_PATHWAYS_CACHE[cache_key] = result
    return result

class CustomHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def do_GET(self):
        if self.path.startswith('/api/summarize'):
            # Parse query parameters
            query = urllib.parse.urlparse(self.path).query
            params = urllib.parse.parse_qs(query)
            gene = params.get('gene', [''])[0]
            name = params.get('name', [''])[0]
            
            # Get API Key from request headers
            api_key = self.headers.get('X-Gemini-API-Key', '')
            
            self.send_response(200)
            self.send_header('Content-Type', 'application/json; charset=utf-8')
            self.send_header('Access-Control-Allow-Origin', '*')
            self.end_headers()
            
            try:
                result = self.perform_summarize(gene, name, api_key)
                self.wfile.write(json.dumps(result).encode('utf-8'))
            except Exception as e:
                self.wfile.write(json.dumps({"error": str(e)}).encode('utf-8'))
        elif self.path.startswith('/api/pathway'):
            # Parse query parameters
            query = urllib.parse.urlparse(self.path).query
            params = urllib.parse.parse_qs(query)
            pathway = params.get('pathway', [''])[0]
            
            # Get API Key from request headers
            api_key = self.headers.get('X-Gemini-API-Key', '')
            
            self.send_response(200)
            self.send_header('Content-Type', 'application/json; charset=utf-8')
            self.send_header('Access-Control-Allow-Origin', '*')
            self.end_headers()
            
            try:
                result = self.perform_pathway_analysis(pathway, api_key)
                self.wfile.write(json.dumps(result).encode('utf-8'))
            except Exception as e:
                self.wfile.write(json.dumps({"error": str(e)}).encode('utf-8'))
        elif self.path.startswith('/api/gene_assistant'):
            # Parse query parameters
            query = urllib.parse.urlparse(self.path).query
            params = urllib.parse.parse_qs(query)
            q_text = params.get('query', [''])[0]
            
            # Get API Key from request headers
            api_key = self.headers.get('X-Gemini-API-Key', '')
            
            self.send_response(200)
            self.send_header('Content-Type', 'application/json; charset=utf-8')
            self.send_header('Access-Control-Allow-Origin', '*')
            self.end_headers()
            
            try:
                result = self.perform_gene_analysis(q_text, api_key)
                self.wfile.write(json.dumps(result).encode('utf-8'))
            except Exception as e:
                self.wfile.write(json.dumps({"error": str(e)}).encode('utf-8'))
        elif self.path.startswith('/api/kegg_pathways'):
            # Parse query parameters
            query = urllib.parse.urlparse(self.path).query
            params = urllib.parse.parse_qs(query)
            cg_locus = params.get('cg', [''])[0]
            cgl_locus = params.get('cgl', [''])[0]
            
            self.send_response(200)
            self.send_header('Content-Type', 'application/json; charset=utf-8')
            self.send_header('Access-Control-Allow-Origin', '*')
            self.end_headers()
            
            try:
                result = get_gene_pathways_and_go(cg_locus, cgl_locus)
                self.wfile.write(json.dumps(result).encode('utf-8'))
            except Exception as e:
                self.wfile.write(json.dumps({"error": str(e)}).encode('utf-8'))
        else:
            super().do_GET()

    def translate_path(self, path):
        parsed = urllib.parse.urlparse(path)
        path_str = parsed.path
        
        # Route requests starting with /data/ to local data/ folder
        if path_str.startswith('/data/'):
            relative_path = path_str[6:] # strip '/data/'
            return os.path.join(os.getcwd(), 'data', relative_path)
            
        # Route other requests to local web/ folder
        relative_path = path_str.lstrip('/')
        if not relative_path:
            relative_path = 'index.html'
        return os.path.join(os.getcwd(), 'web', relative_path)

    def perform_gene_analysis(self, q_text, api_key):
        if not api_key:
            return {"error": "未提供 Gemini API Key。请在右侧详情面板配置您的 API Key。"}
            
        if "DummyKey" in api_key:
            if "抗逆" in q_text or "stress" in q_text.lower():
                return {
                    "summary": "谷氨酸棒状杆菌在面临热激、渗透压、氧化压力等逆境胁迫时，会通过特定的应激反应机制进行自我保护。其中转录因子 SigH (cg0876/Cgl0809) 和氧化应激调节因子 WhiB4 (cg0350/Cgl0339) 扮演了核心的调控作用，启动下游抗逆基因的表达。",
                    "genes": ["cg0350", "cg0876", "cg0409"]
                }
            else:
                return {
                    "summary": f"针对您查询的基因特征 '{q_text}'，AI 识别到了与之最相关的若干个调控与代谢基因，您可以通过下方列表探索它们各自的网络。",
                    "genes": ["cg0350", "cg0876"]
                }
                
        prompt = f"你是一个专业的微生物学 AI 助手，专门研究谷氨酸棒状杆菌 (Corynebacterium glutamicum) ATCC 13032。\n"
        prompt += f"请深度回答并分析关于基因、功能或调控关系的问题：'{q_text}'。\n\n"
        prompt += "请做以下两件事：\n"
        prompt += "1. 提供一段精炼的学术中文总结，解释与该功能或问题相关的基因特征、生物学通路或调控机制（限 200 字以内，排版美观）。\n"
        prompt += "2. 找出与该功能或问题在 C. glutamicum ATCC 13032 中最相关的核心基因的 locus tags（例如 cg0350, cg0814 等）。\n\n"
        prompt += "请严格以 JSON 格式返回，不要带有任何额外的解释文本或 markdown 代码块标记（如 ```json 等），确保返回内容可直接使用 json.loads() 解析。格式如下：\n"
        prompt += '{\n  "summary": "分析与回答内容...",\n  "genes": ["cg0350", "cg0814"]\n}'
        
        models_to_try = ["gemini-2.5-flash", "gemini-1.5-flash", "gemini-2.0-flash"]
        last_err = None
        
        for model_name in models_to_try:
            try:
                gemini_url = f"https://generativelanguage.googleapis.com/v1beta/models/{model_name}:generateContent?key={api_key}"
                payload = {
                    "contents": [{
                        "parts": [{
                            "text": prompt
                        }]
                    }]
                }
                post_data = json.dumps(payload).encode('utf-8')
                gemini_req = urllib.request.Request(
                    gemini_url,
                    data=post_data,
                    headers={'Content-Type': 'application/json'},
                    method='POST'
                )
                
                with urllib.request.urlopen(gemini_req) as gemini_resp:
                    gemini_data = json.loads(gemini_resp.read().decode('utf-8'))
                    text = gemini_data['candidates'][0]['content']['parts'][0]['text'].strip()
                    
                    if text.startswith("```"):
                        text = re.sub(r'^```(?:json)?\s*', '', text)
                        text = re.sub(r'\s*```$', '', text)
                    
                    parsed = json.loads(text)
                    return {
                        "summary": parsed.get("summary", ""),
                        "genes": parsed.get("genes", [])
                    }
            except urllib.error.HTTPError as he:
                try:
                    err_body = he.read().decode('utf-8')
                    err_json = json.loads(err_body)
                    last_err = err_json.get("error", {}).get("message", err_body)
                except Exception:
                    last_err = f"HTTP Error {he.code}: {he.reason}"
                print(f"Model {model_name} failed: {last_err}")
            except Exception as e:
                last_err = str(e)
                print(f"Model {model_name} failed: {last_err}")
                
        raise Exception(f"所有候选模型生成均失败。最后错误: {last_err}")

    def perform_pathway_analysis(self, pathway, api_key):
        if not api_key:
            return {"error": "未提供 Gemini API Key。请在右侧详情面板配置您的 API Key。"}
            
        if "DummyKey" in api_key:
            if "biotin" in pathway.lower() or "生物素" in pathway:
                return {
                    "summary": "生物素（Biotin，维生素 H）合成通路在谷氨酸棒状杆菌中由 bioBFDA 操纵子等基因编码，是参与羧化酶反应的重要辅因子。该通路的调控由生物素蛋白连接酶 BirA 以及合成酶 BioA/BioB 催化。",
                    "genes": ["cg0814", "cg0815", "cg0817"]
                }
            else:
                return {
                    "summary": f"这是一个关于 '{pathway}' 通路的模拟分析总结，识别到相关的调节因子与代谢基因。",
                    "genes": ["cg0350", "cg0409"]
                }
            
        prompt = f"你是一个专业的微生物学 AI 助手，专门研究谷氨酸棒状杆菌 (Corynebacterium glutamicum) ATCC 13032。\n"
        prompt += f"请深度分析代谢通路或生理调控网络：'{pathway}'。\n\n"
        prompt += "请做以下两件事：\n"
        prompt += "1. 提供一段精炼的学术中文总结，描述该通路的生物化学逻辑、关键限速步骤和生理意义（限 200 字以内，排版美观）。\n"
        prompt += "2. 找出该通路在 C. glutamicum ATCC 13032 中关键的所有关联基因的 locus tags（例如 cg0350, cg0814 等）。\n\n"
        prompt += "请严格以 JSON 格式返回，不要带有任何额外的解释文本或 markdown 代码块标记（如 ```json 等），确保返回内容可直接使用 json.loads() 解析。格式如下：\n"
        prompt += '{\n  "summary": "通路的精炼总结...",\n  "genes": ["cg0350", "cg0814"]\n}'
        
        models_to_try = ["gemini-2.5-flash", "gemini-1.5-flash", "gemini-2.0-flash"]
        last_err = None
        
        for model_name in models_to_try:
            try:
                gemini_url = f"https://generativelanguage.googleapis.com/v1beta/models/{model_name}:generateContent?key={api_key}"
                payload = {
                    "contents": [{
                        "parts": [{
                            "text": prompt
                        }]
                    }]
                }
                post_data = json.dumps(payload).encode('utf-8')
                gemini_req = urllib.request.Request(
                    gemini_url,
                    data=post_data,
                    headers={'Content-Type': 'application/json'},
                    method='POST'
                )
                
                with urllib.request.urlopen(gemini_req) as gemini_resp:
                    gemini_data = json.loads(gemini_resp.read().decode('utf-8'))
                    text = gemini_data['candidates'][0]['content']['parts'][0]['text'].strip()
                    
                    # Strip markdown block wrappers if model output included them
                    if text.startswith("```"):
                        text = re.sub(r'^```(?:json)?\s*', '', text)
                        text = re.sub(r'\s*```$', '', text)
                    
                    parsed = json.loads(text)
                    return {
                        "summary": parsed.get("summary", ""),
                        "genes": parsed.get("genes", [])
                    }
            except urllib.error.HTTPError as he:
                try:
                    err_body = he.read().decode('utf-8')
                    err_json = json.loads(err_body)
                    last_err = err_json.get("error", {}).get("message", err_body)
                except Exception:
                    last_err = f"HTTP Error {he.code}: {he.reason}"
                print(f"Model {model_name} failed: {last_err}")
            except Exception as e:
                last_err = str(e)
                print(f"Model {model_name} failed: {last_err}")
                
        raise Exception(f"所有候选模型生成均失败。最后错误: {last_err}")

    def perform_summarize(self, gene, name, api_key):
        # 1. Search PubMed
        term = f'"Corynebacterium glutamicum" AND ({gene}'
        if name and name != "--" and name != gene:
            term += f' OR {name}'
        term += ')'
        
        search_params = {
            "db": "pubmed",
            "term": term,
            "retmode": "json",
            "retmax": 3
        }
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(search_params)
        
        id_list = []
        try:
            req = urllib.request.Request(search_url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req) as resp:
                data = json.loads(resp.read().decode('utf-8'))
                id_list = data.get("esearchresult", {}).get("idlist", [])
        except Exception as e:
            print("PubMed Search Error:", e)
            
        papers = []
        if id_list:
            try:
                fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={','.join(id_list)}&retmode=xml"
                fetch_req = urllib.request.Request(fetch_url, headers={'User-Agent': 'Mozilla/5.0'})
                with urllib.request.urlopen(fetch_req) as fetch_resp:
                    xml_data = fetch_resp.read().decode('utf-8')
                    articles = re.findall(r'<PubmedArticle>(.*?)</PubmedArticle>', xml_data, re.DOTALL)
                    for art in articles:
                        title_match = re.search(r'<ArticleTitle>(.*?)</ArticleTitle>', art, re.DOTALL)
                        abstract_parts = re.findall(r'<AbstractText[^>]*>(.*?)</AbstractText>', art, re.DOTALL)
                        pmid_match = re.search(r'<PMID[^>]*>(.*?)</PMID>', art)
                        
                        title = title_match.group(1).strip() if title_match else "No Title"
                        title = re.sub(r'<[^>]*>', '', title)
                        abstract = " ".join([re.sub(r'<[^>]*>', '', part.strip()) for part in abstract_parts])
                        pmid = pmid_match.group(1).strip() if pmid_match else ""
                        
                        papers.append({
                            "pmid": pmid,
                            "title": title,
                            "abstract": abstract
                        })
            except Exception as e:
                print("PubMed Fetch Error:", e)
                
        # 2. Call Gemini API
        summary = ""
        if not api_key:
            summary = "未提供 Gemini API Key。请在右侧详情面板配置您的 API Key 以生成 AI 智能文献总结。"
        else:
            models_to_try = ["gemini-2.5-flash", "gemini-1.5-flash", "gemini-2.0-flash"]
            last_err = None
            for model_name in models_to_try:
                try:
                    gemini_url = f"https://generativelanguage.googleapis.com/v1beta/models/{model_name}:generateContent?key={api_key}"
                    
                    # Formulate prompt
                    prompt = f"你是一个专业的微生物学 AI 助手，专门研究谷氨酸棒状杆菌 (Corynebacterium glutamicum)。\n"
                    prompt += f"请为基因 {gene} (显示名/常用名: {name if name and name != '--' else '无'}) 生成一份文献与功能总结。\n\n"
                    
                    if papers:
                        prompt += "以下是我们在 PubMed 数据库中检索到的关于该基因的相关研究文献摘要：\n"
                        for idx, paper in enumerate(papers):
                            prompt += f"文献 {idx+1}: {paper['title']}\nPMID: {paper['pmid']}\n摘要: {paper['abstract']}\n\n"
                        prompt += "请根据上述文献的摘要，总结该基因的核心功能、调控机制以及在代谢工程/工业生产中的应用。如果文献中没有涉及某些方面，请结合你所掌握的学术知识进行合理的补充与推断。\n"
                    else:
                        prompt += "我们在 PubMed 中未检索到与该基因直接对应的专属文献。请结合你所掌握的 C. glutamicum 学术知识，详细阐述该基因/转录因子/小RNA 的预测功能、调控通路、以及相关生物学特性。\n"
                    
                    prompt += "\n总结要求：\n1. 使用条理清晰的中文，按以下结构分段总结：【基因概览】、【文献核心研究】、【调控网络与功能】、【发酵应用/科研价值】。\n2. 语言学术、严谨、排版美观（使用 Markdown 格式展示标题和列表）。"
                    
                    payload = {
                        "contents": [{
                            "parts": [{
                                "text": prompt
                            }]
                        }]
                    }
                    
                    post_data = json.dumps(payload).encode('utf-8')
                    gemini_req = urllib.request.Request(
                        gemini_url,
                        data=post_data,
                        headers={'Content-Type': 'application/json'},
                        method='POST'
                    )
                    
                    with urllib.request.urlopen(gemini_req) as gemini_resp:
                        gemini_data = json.loads(gemini_resp.read().decode('utf-8'))
                        summary = gemini_data['candidates'][0]['content']['parts'][0]['text']
                        break
                except urllib.error.HTTPError as he:
                    try:
                        err_body = he.read().decode('utf-8')
                        err_json = json.loads(err_body)
                        last_err = err_json.get("error", {}).get("message", err_body)
                    except Exception:
                        last_err = f"HTTP Error {he.code}: {he.reason}"
                    print(f"Model {model_name} failed: {last_err}")
                except Exception as e:
                    last_err = str(e)
                    print(f"Model {model_name} failed: {last_err}")
            
            if not summary:
                summary = f"Gemini API 总结生成失败。错误信息: {last_err}。\n我们已为您抓取到了相关文献元数据，请参考底部的文献列表。"
                
        return {
            "gene": gene,
            "name": name,
            "summary": summary,
            "papers": [{"pmid": p["pmid"], "title": p["title"]} for p in papers]
        }

class ThreadingHTTPServer(socketserver.ThreadingMixIn, http.server.HTTPServer):
    pass

def open_browser():
    # Wait 1 second to make sure the server has started
    time.sleep(1.0)
    url = f"http://localhost:{PORT}/index.html"
    print(f"Opening network explorer at: {url}")
    webbrowser.open(url)

if __name__ == "__main__":
    server_address = ("", PORT)
    
    # Configure server to allow port re-use (avoid "address already in use" errors)
    socketserver.TCPServer.allow_reuse_address = True
    
    try:
        server = ThreadingHTTPServer(server_address, CustomHTTPRequestHandler)
        print(f"Local Server successfully started on port {PORT}")
        print("Press Ctrl+C to stop the server.")
        
        # Start browser in a background thread
        browser_thread = threading.Thread(target=open_browser)
        browser_thread.daemon = True
        browser_thread.start()
        
        # Serve requests
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping local server. Goodbye!")
        sys.exit(0)
    except Exception as e:
        print(f"Error starting server: {e}")
        sys.exit(1)
